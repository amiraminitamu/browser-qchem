// Scientific regression checks for internal coordinates, figure data and picking policy.
import assert from 'node:assert/strict';
import fs from 'node:fs';
import {parseGeometry,WATER_XYZ,ANGSTROM_PER_BOHR,toXYZ} from '../web/molecule.mjs';
import {coordinateValue,coordinateLabel} from '../web/coordinates.mjs';
import {scanSpec,scanGeometry,distance,runScan,scanCSV,scanXYZ} from '../web/scan.mjs';
import {pesFigureModel,figureSVG,imagePDF} from '../web/pes-export.mjs';
import {MoleculeViewer} from '../web/viewer.mjs';
import {TASKS} from '../web/workbench.mjs';
const near=(a,b,t=1e-9)=>assert.ok(Math.abs(a-b)<t,`${a} != ${b}`);
const water=parseGeometry(WATER_XYZ),waterBefore=structuredClone(water);
const bend=scanSpec(water,{coordinate:'angle',atoms:[1,0,2],moving:[2],start:80,end:180,count:5,unit:'degree'});
for(const target of [80,104.52,140,180]){
  const geometry=scanGeometry(water,bend,target);
  near(coordinateValue(geometry,'angle',[1,0,2]),target,1e-6);
  near(distance(geometry,0,2),distance(water,0,2));
  assert.deepEqual(geometry.atoms[0],water.atoms[0]);assert.deepEqual(geometry.atoms[1],water.atoms[1]);
}
assert.deepEqual(water,waterBefore);
// Bend a rigid multi-atom fragment without changing any of its internal distances.
const fragment=parseGeometry('H 0 1 0\nC 0 0 0\nC 1.5 0 0\nH 2.2 .8 .3');
const fspec=scanSpec(fragment,{coordinate:'angle',atoms:[0,1,2],moving:[2,3],start:30,end:170,count:3});
for(const target of [30,90,170]){const g=scanGeometry(fragment,fspec,target);near(coordinateValue(g,'angle',fspec.atoms),target);near(distance(g,2,3),distance(fragment,2,3));}
const peroxide=parseGeometry('H -.3 .9 0\nO 0 0 0\nO 1.45 0 0\nH 1.75 0 .9');
const torsion=scanSpec(peroxide,{coordinate:'dihedral',atoms:[0,1,2,3],moving:[2,3],start:0,end:360,count:7});
const wrap=angle=>((angle+180)%360+360)%360-180;
for(const target of [-180,-90,0,30,180,270,360,720]){
  const g=scanGeometry(peroxide,torsion,target);
  near(wrap(coordinateValue(g,'dihedral',torsion.atoms)-target),0);
  for(const [a,b] of [[0,1],[1,2],[2,3]])near(distance(g,a,b),distance(peroxide,a,b));
  near(coordinateValue(g,'angle',[0,1,2]),coordinateValue(peroxide,'angle',[0,1,2]));
  near(coordinateValue(g,'angle',[1,2,3]),coordinateValue(peroxide,'angle',[1,2,3]));
  assert.deepEqual(g.atoms[0],peroxide.atoms[0]);assert.deepEqual(g.atoms[1],peroxide.atoms[1]);
}
const end=scanGeometry(peroxide,torsion,360),start=scanGeometry(peroxide,torsion,0);
end.atoms.forEach((atom,i)=>atom.xyz.forEach((x,d)=>near(x,start.atoms[i].xyz[d])));
const linear=parseGeometry('H -1 0 0\nO 0 0 0\nH 1 0 0');
const linearSpec=scanSpec(linear,{coordinate:'angle',atoms:[0,1,2],moving:[2],start:150,end:180,count:3});
near(coordinateValue(scanGeometry(linear,linearSpec,150),'angle',linearSpec.atoms),150);
assert.throws(()=>coordinateValue(linear,'dihedral',[0,1,2,0]),/collinear/);
for(const bad of [{atoms:[1,0,1]},{start:0},{end:181},{moving:[0,2]},{unit:'angstrom'}])assert.throws(()=>scanSpec(water,{...bend,...bad}));
for(const bad of [{coordinate:'other'},{atoms:[0,1,2,2]},{moving:[0,3]},{unit:'bohr'}])assert.throws(()=>scanSpec(peroxide,{...torsion,...bad}));
assert.match(coordinateLabel(water,bend),/θ\(H₂–O₁–H₃\) \[°\]/);
console.log('Bends, signed torsions, periodicity, fixed atoms, rigid fragments, singularities and validation PASS');

// Exercise actual WASM for both new coordinate types, with independent references.
const createQChem=(await import('../build/wasm/dist/qc_core.mjs')).default,engine=await createQChem();
const angleScan=await runScan(engine,water,{...bend,start:90,end:120,count:3},{method:'rhf'});
assert.equal(angleScan.completed,3);assert.ok(angleScan.points.every(p=>p.converged));
for(const p of angleScan.points){near(coordinateValue(p.molecule,'angle',bend.atoms),p.coordinate);assert.equal(p.coordinateBohr,undefined);assert.equal(p.coordinateDegrees,p.coordinate);}
const rotation=await runScan(engine,peroxide,{...torsion,start:30,end:150,count:3},{method:'rhf'});
assert.ok(rotation.points.every(p=>p.converged));
for(const p of rotation.points)near(wrap(coordinateValue(p.molecule,'dihedral',torsion.atoms)-p.coordinate),0);
const raw=scanCSV(rotation);assert.match(raw,/coordinate_type,atom_indices/);assert.match(raw,/"dihedral","1 2 3 4"/);
const xyz=scanXYZ(rotation).trim().split('\n');assert.equal(xyz.length,18);
for(let i=0;i<3;i++)near(wrap(coordinateValue(parseGeometry(xyz.slice(i*6,i*6+6).join('\n')),'dihedral',torsion.atoms)-rotation.points[i].coordinate),0);
console.log('Actual-WASM water bend and peroxide torsion scans; degree-aware CSV and XYZ round trips PASS');

const saved={...angleScan,source:water,moleculeName:'Water',formula:'H₂O'};
const model=pesFigureModel(saved,{unit:'ev',relative:true});
assert.equal(Math.min(...model.points.map(p=>p.y)),0);assert.match(model.xLabel,/θ/);assert.match(model.yLabel,/eV/);
const absolute=pesFigureModel(saved,{unit:'eh',relative:false});near(absolute.points[0].y,saved.points[0].energy);
const kcal=pesFigureModel(saved,{unit:'kcal'});near(kcal.points[0].y,model.points[0].y*627.5094740631/27.211386245988);
const withFailure=structuredClone(saved);withFailure.points[0].converged=false;withFailure.points[0].energy=-1000;
assert.notEqual(pesFigureModel(withFailure).minimum.index,0);
withFailure.points[0].energy=null;assert.equal(pesFigureModel(withFailure).points[0].y,null);
const malicious=pesFigureModel({...saved,moleculeName:'<script>&"'});assert.ok(!figureSVG(malicious).svg.includes('<script>'));
const figure=figureSVG(model);assert.match(figure.svg,/fill="#ffffff"/);assert.match(figure.svg,/SAMPLED MINIMUM/);assert.match(figure.svg,/H₂–O₁–H₃/);
const pdf=imagePDF(new Uint8Array([255,216,255,217]),2700,2100),text=new TextDecoder().decode(pdf);
assert.ok(text.startsWith('%PDF-1.4'));const xref=Number(text.match(/startxref\n(\d+)/)[1]);assert.equal(text.slice(xref,xref+4),'xref');
assert.match(text,/\/Width 2700 \/Height 2100/);
if(process.env.PES_QA_DIR){fs.mkdirSync(process.env.PES_QA_DIR,{recursive:true});fs.writeFileSync(process.env.PES_QA_DIR+'/figure.svg',figure.svg);fs.writeFileSync(process.env.PES_QA_DIR+'/model.json',JSON.stringify(model));}
console.log('Energy conversion, converged minimum, missing energies, caption escaping, metadata and PDF byte offsets PASS');

// The adapter must remove all picking affordances when leaving PES.
const calls={clickable:[],hoverable:[],shapes:[],labels:[]};let selected=0;
const api={setProjection(){},resize(){},removeAllModels(){},getView:()=>[],setView(){},rotate(){},zoomTo(){},zoom(){},modelToScreen:()=>[],render(){},setStyle(){},
  addModel:()=>({addAtoms(){},selectedAtoms:()=>[]}),setClickable:(_s,on,fn)=>calls.clickable.push({on,fn}),setHoverable:(_s,on)=>calls.hoverable.push(on),
  removeShape:shape=>{calls.shapes=calls.shapes.filter(s=>s!==shape);},removeLabel:label=>{calls.labels=calls.labels.filter(l=>l!==label);},
  addLabel:(text,spec)=>{const label={text,spec};calls.labels.push(label);return label;}};
for(const method of ['addLine','addSphere','addCylinder','addCustom','addArrow'])api[method]=spec=>{const shape={method,spec,updateStyle(){}};calls.shapes.push(shape);return shape;};
globalThis.$3Dmol={createViewer:()=>api};globalThis.ResizeObserver=class{observe(){}};
const viewer=new MoleculeViewer({dataset:{},clientWidth:800,clientHeight:700},()=>selected++);
viewer.setMolecule(water);viewer.select([0,1]);
assert.equal(calls.clickable.at(-1).on,false);assert.equal(calls.hoverable.at(-1),false);assert.equal(calls.shapes.length,0);
viewer.setPickingMode('internal-coordinate','angle');viewer.select([1,0,2]);assert.equal(calls.clickable.at(-1).on,true);assert.ok(calls.shapes.some(s=>s.method==='addCustom'));
const oldCallback=calls.clickable.at(-1).fn;
viewer.setPickingMode('none');oldCallback({index:1});assert.equal(selected,0);assert.equal(calls.shapes.length,0);assert.equal(calls.labels.some(l=>/^\d+:/.test(l.text)),false);
for(const task of TASKS)assert.equal(task.picking,task.id==='scan'?'internal-coordinate':'none');
console.log('Clean inspection modes disable callbacks, hover and gizmos, including stale picking callbacks PASS');
