// Actual WASM tests for the basis catalog, scan geometry, cancellation and exports.
import assert from 'node:assert/strict';
import fs from 'node:fs';
import createQChem from '../build/wasm/dist/qc_core.mjs';
import {parseGeometry,packBasis,basisSummary,toXYZ,ANGSTROM_PER_BOHR} from '../web/molecule.mjs';
import {computeRHF as originalComputeRHF} from '../web/rhf.mjs';
import {scanSpec,scanGeometry,runScan as originalRunScan,distance,connectedFragment,scanCSV,scanXYZ} from '../web/scan.mjs';
const engine=await createQChem();
const reference=JSON.parse(fs.readFileSync(new URL('./fixtures/workspace.json',import.meta.url)));
const near=(actual,expected,tol=1e-9)=>assert.ok(Math.abs(actual-expected)<tol,`${actual} != ${expected}`);
for(const basis of ['sto-3g','6-31g','cc-pvdz']) {
  const water=reference.cases.find(c=>c.name==='water'&&c.basis===basis);
  const molecule=parseGeometry(water.geometry,water.unit,basis);
  const result=await computeRHF(engine,molecule);
  assert.equal(result.nao,water.nao);assert.ok(result.rhf.converged);near(result.rhf.energy,water.energy);
  const h2=parseGeometry('H 0 0 0\nH 0 0 .74','angstrom',basis);
  const scan=await runScan(engine,h2,{a:0,b:1,moving:[1],start:.5,end:1.5,count:5});
  assert.equal(scan.completed,5);assert.equal(scan.cancelled,false);
  for(const point of scan.points) {
    assert.ok(point.converged);assert.equal(point.result.eri,undefined);
    near(distance(point.molecule,0,1)*ANGSTROM_PER_BOHR,point.coordinate,1e-12);
    const ref=reference.cases.find(c=>c.basis===basis&&c.name.startsWith('h2-')&&Number(c.name.slice(3))===point.coordinate);
    near(point.energy,ref.energy);
    assert.deepEqual(point.molecule.atoms[0].xyz,h2.atoms[0].xyz);
  }
  const csv=scanCSV(scan).trim().split('\n');assert.equal(csv.length,6);assert.ok(csv[0].includes('energy_hartree'));
  const xyz=scanXYZ(scan).trim().split('\n');assert.equal(xyz.length,20);
  for(let i=0;i<5;i++)near(distance(parseGeometry(xyz.slice(4*i,4*i+4).join('\n')),0,1)*ANGSTROM_PER_BOHR,.5+.25*i,1e-11);
  console.log(`${basis}: water and five-point H2 PES match PySCF; CSV/XYZ exports PASS`);
}
const fragment=parseGeometry('C 0 0 0\nC 1.54 0 0\nH -1 0 0\nH 2.54 0 0');
assert.deepEqual(connectedFragment(fragment,0,1),[1,3]);
const spec=scanSpec(fragment,{a:0,b:1,moving:[1,3],start:4,end:2,count:3,unit:'bohr'});
for(const target of [4,3,2]) {
  const geometry=scanGeometry(fragment,spec,target);near(distance(geometry,0,1),target,1e-12);
  near(distance(geometry,1,3),distance(fragment,1,3),1e-12);
  assert.deepEqual(geometry.atoms[0],fragment.atoms[0]);assert.deepEqual(geometry.atoms[2],fragment.atoms[2]);
}
assert.deepEqual(parseGeometry(toXYZ(fragment)).atoms.map(a=>a.symbol),fragment.atoms.map(a=>a.symbol));
const ring=parseGeometry('C 0 0 0\nC 1.4 0 0\nC .7 1.212 0');
assert.throws(()=>connectedFragment(ring,0,1),/ring/);
for(const bad of [{a:0,b:0},{moving:[0,1]},{moving:[2]},{count:102},{start:0},{start:1,end:1},{unit:'nm'}])
  assert.throws(()=>scanSpec(fragment,{a:0,b:1,start:1,end:2,count:3,...bad}));
const tooBig=parseGeometry(Array.from({length:32},(_,i)=>`Ne 0 0 ${i*3}`).join('\n'),'angstrom','cc-pvdz');
assert.equal(basisSummary(tooBig).fits,false);assert.throws(()=>packBasis(tooBig),/256 AOs/);
assert.throws(()=>parseGeometry('H 0 0 0','angstrom','fake'),/basis/);
const h2=parseGeometry('H 0 0 0\nH 0 0 .74');
const config={a:0,b:1,start:.5,end:1.5,count:5};
let cancelled=false;
const stopped=await runScan(engine,h2,config,{isCancelled:()=>cancelled},()=>cancelled=true);
assert.equal(stopped.completed,1);assert.equal(stopped.cancelled,true);assert.equal(scanCSV(stopped).trim().split('\n').length,2);
const empty=await runScan(engine,h2,config,{isCancelled:()=>true});assert.equal(empty.completed,0);assert.equal(empty.cancelled,true);
let iterations=0;
const midSCF=await runScan(engine,h2,config,{isCancelled:()=>iterations>0},()=>{},()=>iterations++);
assert.equal(midSCF.completed,0);assert.equal(midSCF.cancelled,true);
const water=parseGeometry(reference.cases.find(c=>c.name==='water').geometry);
const limited=await runScan(engine,water,config,{maxIterations:1});assert.ok(limited.points.every(p=>p.status==='not_converged'&&!p.converged));
const collision=parseGeometry('H 0 0 0\nH 1 0 0\nH 2 0 0');
const failed=await runScan(engine,collision,{a:0,b:1,start:1,end:3,count:3},{charge:1});
assert.equal(failed.completed,3);assert.equal(failed.points[1].status,'error');assert.equal(failed.points[1].energy,null);assert.match(failed.points[1].error,/coincide/);
assert.ok(failed.points[2].result);assert.equal(failed.points[2].status,'converged');
await assert.rejects(()=>runScan(engine,h2,config,{charge:1}),/even electron/);
const recovered=await computeRHF(engine,h2);assert.ok(recovered.rhf.converged);
console.log('Rigid fragments, reverse ranges, bohr units, rings, AO limits, cancellation, failed points, recovery PASS');

function computeRHF(e,m,o={},...rest){return originalComputeRHF(e,m,{...o,scfIntegrals:'exact'},...rest);}

function runScan(e,m,c,o={},...rest){return originalRunScan(e,m,c,{...o,scfIntegrals:'exact'},...rest);}
