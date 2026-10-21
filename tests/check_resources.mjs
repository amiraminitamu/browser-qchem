// SPDX-License-Identifier: AGPL-3.0-only
import assert from 'node:assert/strict';
import fs from 'node:fs';
import {parseGeometry,packBasis,basisSummary,WATER_XYZ} from '../web/molecule.mjs';
import {estimateResources,assertResources,deviceBudget,jobBudget,retainedBytes} from '../web/resources.mjs';
import {MIB,MAX_AO} from '../web/limits.mjs';
import {computeElectronic,validateElectronicOptions} from '../web/electronic.mjs';
import {computeIntegrals,computeERI} from '../web/engine.mjs';
import {runScan} from '../web/scan.mjs';
import createQChem from '../build/wasm/dist/qc_core.mjs';
const water=parseGeometry(WATER_XYZ),chain=n=>parseGeometry(Array.from({length:n},(_,i)=>`H 0 0 ${i*.85}`).join('\n'),'angstrom','cc-pvdz');
assert.equal(MAX_AO,256);assert.equal(basisSummary(chain(20)).nao,100);assert.equal(packBasis(chain(20)).aoLabels.length,100);
assert.equal(deviceBudget({}),256*MIB);assert.equal(deviceBudget({deviceMemory:8}),512*MIB);assert.equal(deviceBudget({deviceMemory:2}),128*MIB);
for(const v of [NaN,-1,Infinity,undefined])assert.equal(jobBudget(v),256*MIB);
assert.equal(jobBudget(999999*MIB),512*MIB);assert.equal(jobBudget(1),64*MIB);
const normal=estimateResources(water);assert.equal(normal.blocked,false);assert.equal(normal.level,'light');
const larger=estimateResources(chain(20));assert.equal(larger.blocked,false);assert.equal(larger.nao,100);assert.equal(larger.level,'heavy');
const df=estimateResources(chain(12)),exact=estimateResources(chain(12),{scfIntegrals:'exact'}),mp2=estimateResources(chain(12),{method:'mp2'}),ri=estimateResources(chain(12),{method:'df-rmp2'});
assert.ok(exact.peakBytes>df.peakBytes);assert.ok(ri.correlationFunctions>0);assert.ok(mp2.peakBytes>=df.peakBytes);assert.ok(mp2.workScore>df.workScore);
assert.equal(estimateResources(chain(20),{scfIntegrals:'exact'}).blocked,true);
assert.equal(estimateResources(chain(20),{mode:'integrals'}).blocked,true);
assert.equal(estimateResources(chain(28)).blocked,true);
assert.equal(estimateResources(chain(28),{resourceBudgetBytes:512*MIB}).blocked,false);
const few=estimateResources(chain(20),{mode:'scan',scanCount:3}),many=estimateResources(chain(20),{mode:'scan',scanCount:101});
assert.ok(many.peakBytes>few.peakBytes);assert.ok(many.workScore>few.workScore);assert.equal(many.blocked,true);
assert.match(many.summary,/101 sequential points/);assert.ok(many.peakBytes<101*larger.peakBytes,'Sequential scans must not multiply the full per-point workspace');
const manyElectrons=parseGeometry(Array.from({length:10},(_,i)=>`Ne 0 0 ${i*4}`).join('\n'));
assert.equal(validateElectronicOptions(manyElectrons).electrons,100);
assert.equal(validateElectronicOptions(manyElectrons,{method:'uhf'}).nalpha,50);
const tooLarge=parseGeometry(Array.from({length:32},(_,i)=>`Ne 0 0 ${i*3}`).join('\n'),'angstrom','cc-pvdz');
assert.equal(basisSummary(tooLarge).fits,false);assert.throws(()=>packBasis(tooLarge),/256 AOs/);
assert.throws(()=>assertResources(tooLarge),/256 AOs/);assert.throws(()=>validateElectronicOptions(chain(20),{scfIntegrals:'exact'}),/DF-JK/);
const buffer=new Float64Array(1024);assert.equal(retainedBytes({a:buffer,b:buffer.subarray(1),c:{a:buffer}}),8192);
console.log('Backend-aware peak memory, auxiliary/work estimates, sequential scan retention, device budgets, structural/native ceilings and unique history buffers PASS');
const engine=await createQChem(),fixture=JSON.parse(fs.readFileSync(new URL('./fixtures/resource-large.json',import.meta.url)));
const near=(a,b,t=3e-8)=>assert.ok(Math.abs(a-b)<t,`${a} != ${b}`);
for(const c of fixture.cases){
  const molecule=parseGeometry(c.geometry,c.unit,c.basis),r=await computeElectronic(engine,molecule,c.options),ref=r.uhf??r.rhf;
  assert.ok(r.nao>40);assert.equal(ref.converged,true);assert.equal(r.eri,undefined);near(ref.energy,c.referenceEnergy,3e-9);
  for(const key of ['energy','correlationEnergy','oppositeSpinEnergy','sameSpinEnergy'])near(r.mp2[key],c[key]);
  assert.equal(r.scfIntegrals.auxiliaryFunctions,c.jkfitFunctions);near(computeERI(engine,molecule,c.eriIndices),c.exactERI,2e-10);
  if(c.conventionalCorrelationEnergy!==null){const conventional=await computeElectronic(engine,molecule,{...c.options,method:'mp2'});near(conventional.mp2.correlationEnergy,c.conventionalCorrelationEnergy);assert.equal(conventional.eri,undefined);}
}
// Exact remains available above 40 when its own resource estimate fits.
const e=fixture.exactCase,neon=parseGeometry(e.geometry,e.unit,e.basis);
const exactResult=await computeElectronic(engine,neon,{scfIntegrals:'exact'});near(exactResult.rhf.energy,e.energy,2e-9);assert.equal(exactResult.nao,42);assert.ok(exactResult.rhf.converged);assert.equal(exactResult.eri.length,42**4);
let allocations=0;const allocating=engine._malloc;engine._malloc=(...args)=>{++allocations;return allocating(...args);};
await assert.rejects(()=>computeElectronic(engine,chain(20),{scfIntegrals:'exact'}),/exceeds|budget/);
assert.throws(()=>computeIntegrals(engine,chain(20)),/exceeds|budget/);
await assert.rejects(()=>runScan(engine,chain(20),{a:0,b:1,start:.8,end:1,count:101}),/budget/);
assert.equal(allocations,0,'Unsafe jobs must fail before any WASM input or integral allocation');engine._malloc=allocating;
console.log('Real-WASM 50/72 AO RHF/UHF, both DF-MP2 flavors, conventional MP2, exact 42 AO SCF, on-demand ERIs and allocation-free over-budget rejection PASS');

assert.equal(engine._malloc(600*MIB),0,'A request beyond the configured WASM heap must fail without aborting');
assert.ok((await computeElectronic(engine,water)).rhf.converged);
console.log('512 MiB WASM heap ceiling returns allocation failure without aborting; subsequent calculation recovers PASS');
