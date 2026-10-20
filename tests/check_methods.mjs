// Actual WASM numerical regressions, cancellation and method-aware scans.
import assert from 'node:assert/strict';
import fs from 'node:fs';
import createQChem from '../build/wasm/dist/qc_core.mjs';
import {parseGeometry,WATER_XYZ} from '../web/molecule.mjs';
import {computeElectronic,totalEnergy} from '../web/electronic.mjs';
import {serializable} from '../web/engine.mjs';
import {runScan,scanCSV,scanXYZ} from '../web/scan.mjs';
const engine=await createQChem(),water=parseGeometry(WATER_XYZ),h2=parseGeometry('H 0 0 0\nH 0 0 3');
const near=(a,b,tol=2e-8)=>assert.ok(Math.abs(a-b)<tol,`${a} != ${b}`);
const fixture=JSON.parse(fs.readFileSync(new URL('./fixtures/methods.json',import.meta.url)));
for(const c of fixture.cases) {
  const result=await computeElectronic(engine,parseGeometry(c.geometry,c.unit,c.basis),c.options),r=result.uhf??result.rhf;
  assert.ok(r.converged);near(r.energy,c.referenceEnergy);
  if(result.uhf)near(r.spinSquare,c.spinSquare,2e-6);
  else for(const key of ['energy','correlationEnergy','oppositeSpinEnergy','sameSpinEnergy'])near(result.mp2[key],c[key]);
  const output=serializable(result);assert.ok(Array.isArray(output.overlap));
  if(result.uhf)assert.ok(Array.isArray(output.uhf.densityAlpha));
}
console.log(`${fixture.cases.length} UHF/MP2 PySCF fixtures and JSON arrays PASS`);
const rhf=await computeElectronic(engine,h2),broken=await computeElectronic(engine,h2,{method:'uhf'}),core=await computeElectronic(engine,h2,{method:'uhf',guess:'core'});
assert.ok(broken.uhf.energy<rhf.rhf.energy-.1);assert.ok(broken.uhf.spinSquare>.9);near(core.uhf.energy,rhf.rhf.energy);assert.equal(core.uhf.mixedGuessApplied,false);
const noVirtual=await computeElectronic(engine,parseGeometry('He 0 0 0'),{method:'mp2'});assert.equal(noVirtual.mp2.correlationEnergy,0);assert.equal(noVirtual.mp2.virtual,0);assert.equal(noVirtual.mp2.minimumAbsoluteDenominator,null);
const limited=await computeElectronic(engine,water,{method:'mp2',maxIterations:1});assert.equal(limited.mp2.status,'not_computed');assert.equal(totalEnergy(limited),null);
for(const options of [{method:'ump2'},{method:'uhf',multiplicity:2},{method:'uhf',multiplicity:0},{method:'uhf',multiplicity:1.1},{method:'uhf',multiplicity:99},{method:'uhf',charge:11},{method:'uhf',charge:.1},{method:'uhf',guess:'bad'},{method:'uhf',maxIterations:0},{method:'mp2',charge:1}])await assert.rejects(()=>computeElectronic(engine,water,options));
await assert.rejects(()=>computeElectronic(engine,parseGeometry('H 0 0 0'),{method:'uhf'}),/Multiplicity/);
await assert.rejects(()=>computeElectronic(engine,parseGeometry('H 0 0 0\nH 0 0 .00001','bohr'),{method:'uhf'}),/linearly dependent/);
assert.ok((await computeElectronic(engine,water,{method:'uhf',guess:'core',diis:false})).uhf.converged);
const ul=await computeElectronic(engine,water,{method:'uhf',maxIterations:1});assert.equal(ul.uhf.converged,false);
let stop=false,blocks=0;
await assert.rejects(()=>computeElectronic(engine,water,{method:'mp2',isCancelled:()=>stop},()=>{},()=>{blocks++;stop=true;}),{name:'AbortError'});assert.equal(blocks,1);
stop=false;await assert.rejects(()=>computeElectronic(engine,water,{method:'uhf',isCancelled:()=>stop},()=>{stop=true;}),{name:'AbortError'});
for(const method of ['uhf','mp2']) {
 const config={a:0,b:1,start:.7,end:1.5,count:3,unit:'angstrom'};
 const scan=await runScan(engine,h2,config,{method,maxIterations:500});assert.equal(scan.completed,3);assert.equal(scan.spec.method,method);
 for(const point of scan.points) {assert.ok(point.converged);near(point.energy,totalEnergy(point.result),1e-12);assert.equal(point.result.eri,undefined);}
 assert.match(scanCSV(scan),/spin_square,correlation_energy_hartree/);assert.match(scanXYZ(scan),new RegExp(`method=${method.toUpperCase()}`));
 let stopped=false;const partial=await runScan(engine,h2,config,{method,isCancelled:()=>stopped},()=>stopped=true);assert.equal(partial.completed,1);assert.ok(partial.cancelled);
}
const failed=await runScan(engine,water,{a:0,b:1,start:.7,end:1.5,count:2},{method:'mp2',maxIterations:1});assert.ok(failed.points.every(p=>p.energy===null&&!p.converged));
assert.equal(engine._qc_mp2_create(0),0);assert.equal(engine._qc_mp2_step(0),-1);assert.equal(engine._qc_uhf_step(0),-1);assert.equal(engine._qc_uhf_create(0,1,0,100,1,0,1e-10,1e-8,1e-8),0);
assert.ok((await computeElectronic(engine,water,{method:'mp2'})).mp2.energy<-74.99);
console.log('Broken symmetry, spin validation, no virtuals, failed references, cancellation, UHF/MP2 scans and recovery PASS');
