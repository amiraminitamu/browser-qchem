// SPDX-License-Identifier: AGPL-3.0-only
// Real WASM against independent PySCF RI-MP2 energies and spin decompositions.
import assert from 'node:assert/strict';
import fs from 'node:fs';
import createQChem from '../build/wasm/dist/qc_core.mjs';
import {parseGeometry,WATER_XYZ} from '../web/molecule.mjs';
import {computeElectronic,totalEnergy,validateElectronicOptions} from '../web/electronic.mjs';
import {runScan,scanCSV,scanXYZ} from '../web/scan.mjs';
import {serializable} from '../web/engine.mjs';
const engine=await createQChem();
const fixture=JSON.parse(fs.readFileSync(new URL('./fixtures/df-mp2.json',import.meta.url)));
const near=(a,b,t=2e-8)=>assert.ok(Math.abs(a-b)<t,`${a} != ${b}`);
for(const c of fixture.cases){
  const m=parseGeometry(c.geometry,c.unit,c.basis),r=await computeElectronic(engine,m,{...c.options,scfIntegrals:'exact'}),p=r.mp2;
  assert.equal(r.method,c.options.method.toUpperCase());assert.equal(p.reference,c.options.method==='df-ump2'?'UHF':'RHF');
  assert.equal(p.status,'complete');assert.equal(p.densityFitted,true);assert.equal(p.allElectron,true);assert.equal(p.auxiliaryBasis,'cc-pVDZ-RI');assert.equal(p.auxiliaryFunctions,c.naux);
  assert.ok(p.metricRank>0&&p.metricRank<=p.auxiliaryFunctions);assert.equal(p.metricEigenvalueCutoff,1e-10);
  near((r.uhf??r.rhf).energy,c.referenceEnergy);
  for(const key of ['energy','correlationEnergy','oppositeSpinEnergy','sameSpinEnergy'])near(p[key],c[key]);
  near(p.sameSpinAlphaEnergy+p.sameSpinBetaEnergy,p.sameSpinEnergy,1e-12);near(p.referenceEnergy+p.correlationEnergy,totalEnergy(r),1e-12);
  assert.equal(serializable(r).mp2.densityFitted,true);
  if((r.uhf??r.rhf).electrons===1){assert.equal(p.correlationEnergy,0);assert.equal(p.minimumAbsoluteDenominator,null);}
}
console.log(`${fixture.cases.length} actual-WASM DF-RMP2/DF-UMP2 PySCF comparisons and spin decomposition PASS`);
const water=parseGeometry(WATER_XYZ),h2=parseGeometry('H 0 0 0\nH 0 0 .74');
const r=await computeElectronic(engine,water,{method:'df-rmp2'}),u=await computeElectronic(engine,water,{method:'df-ump2',guess:'core'});
near(r.mp2.energy,u.mp2.energy);near(r.mp2.oppositeSpinEnergy,u.mp2.oppositeSpinEnergy);near(r.mp2.sameSpinEnergy,u.mp2.sameSpinEnergy);
for(const method of ['df-rmp2','df-ump2']){
  const failed=await computeElectronic(engine,water,{method,maxIterations:1});assert.equal(failed.mp2.status,'not_computed');assert.equal(totalEnergy(failed),null);
  let cancel=false;
  await assert.rejects(()=>computeElectronic(engine,water,{method,isCancelled:()=>cancel},()=>{},p=>{if(p.phase==='metric')cancel=true;}),{name:'AbortError'});
  await assert.rejects(()=>computeElectronic(engine,water,{method,isCancelled:()=>true}),{name:'AbortError'});
  let blocks=0;
  await assert.rejects(()=>computeElectronic(engine,water,{method,isCancelled:()=>blocks>0},()=>{},p=>{if(p.phase==='correlation')++blocks;}),{name:'AbortError'});
  const scan=await runScan(engine,h2,{a:0,b:1,start:.6,end:1,count:3},{method,guess:'core'});
  assert.equal(scan.completed,3);assert.equal(scan.spec.auxiliaryBasis,'cc-pvdz-ri');assert.ok(scan.points.every(p=>p.converged&&p.result.mp2.densityFitted));
  assert.match(scanCSV(scan),/auxiliary_basis,auxiliary_functions,metric_rank,density_fitted/);assert.match(scanCSV(scan),/cc-pVDZ-RI/);assert.match(scanXYZ(scan),/auxiliary_basis=cc-pvdz-ri/);
}
await assert.rejects(()=>computeElectronic(engine,h2,{method:'df-rmp2',charge:1}),/valid RHF reference/);
await assert.rejects(()=>computeElectronic(engine,h2,{method:'df-ump2',charge:1,multiplicity:1}),/Multiplicity/);
assert.throws(()=>validateElectronicOptions(water,{method:'df-rmp2',auxiliaryBasis:'not-a-basis'}),/cc-pVDZ-RI/);
const ion=await runScan(engine,h2,{a:0,b:1,start:.7,end:1.3,count:3},{method:'df-ump2',charge:1,multiplicity:2,guess:'core'});
assert.ok(ion.points.every(p=>p.converged&&p.result.mp2.correlationEnergy===0));assert.match(scanCSV(ion),/DF-UMP2/);
for(const method of ['rhf','uhf','mp2','df-rmp2','df-ump2']){const result=await computeElectronic(engine,water,{method});assert.ok((result.rhf??result.uhf).converged);}
console.log('Restricted/unrestricted limit, failed references, validation, cancellation during fitting and contraction, PES metadata, one-electron ions, repeated recovery PASS');
