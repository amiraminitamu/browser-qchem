// SPDX-License-Identifier: AGPL-3.0-only
// Real WASM compared with independent PySCF DF-JK reference fixtures.
import assert from 'node:assert/strict';
import fs from 'node:fs';
import createQChem from '../build/wasm/dist/qc_core.mjs';
import {parseGeometry,WATER_XYZ} from '../web/molecule.mjs';
import {computeElectronic,validateElectronicOptions} from '../web/electronic.mjs';
import {prepareIntegrals,computeERI} from '../web/engine.mjs';
import {runScan,scanCSV,scanXYZ} from '../web/scan.mjs';
import {pesFigureModel} from '../web/pes-export.mjs';
const engine=await createQChem(),fixture=JSON.parse(fs.readFileSync(new URL('./fixtures/df-scf.json',import.meta.url)));
const near=(a,b,t=2e-8)=>assert.ok(Math.abs(a-b)<t,`${a} != ${b}`);
let dense=0;
const compute=engine._qc_compute_mode;
engine._qc_compute_mode=(flag,...args)=>{dense+=flag;return compute(flag,...args);};
for(const c of fixture.cases){
  const molecule=parseGeometry(c.geometry,c.unit,c.basis),r=await computeElectronic(engine,molecule,c.options),fit=r.scfIntegrals,ref=r.uhf??r.rhf;
  assert.equal(fit.type,'df');assert.equal(fit.densityFitted,true);assert.equal(fit.auxiliaryBasis,'def2-universal-JKFIT');
  assert.equal(fit.auxiliaryFunctions,c.jkfitFunctions);assert.equal(fit.metricRank,c.jkfitRank);assert.equal(fit.metricEigenvalueCutoff,1e-10);
  assert.equal(r.eri,undefined);assert.equal(r.checks.eriSymmetryError,null);assert.ok(r.checks.passed);assert.ok(ref.converged);
  near(ref.energy,c.referenceEnergy,2e-9);
  for(const key of ['energy','correlationEnergy','oppositeSpinEnergy','sameSpinEnergy'])near(r.mp2[key],c[key]);
  if(r.uhf)near(r.uhf.spinSquare,c.spinSquare,2e-7);
  near(computeERI(engine,molecule,c.eriIndices),c.exactERI,2e-10);
  if(c.conventionalCorrelationEnergy!==null){
    const conventional=await computeElectronic(engine,molecule,{...c.options,method:'mp2'});
    assert.equal(conventional.eri,undefined);assert.equal(conventional.mp2.densityFitted,false);
    near(conventional.mp2.correlationEnergy,c.conventionalCorrelationEnergy);
  }
}
assert.equal(dense,0,'DF SCF, exact on-demand inspection, and direct conventional MP2 must not request dense ERIs');
console.log(`${fixture.cases.length} actual-WASM DF-JK + DF-MP2 PySCF comparisons; direct conventional MP2; exact ERI queries; no dense ERI allocation PASS`);
const water=parseGeometry(WATER_XYZ),h2=parseGeometry('H 0 0 0\nH 0 0 .74');
assert.equal((await computeElectronic(engine,water)).scfIntegrals.type,'df');
for(const method of ['rhf','uhf','mp2','df-rmp2','df-ump2']){
  const exact=await computeElectronic(engine,water,{method,scfIntegrals:'exact'});
  assert.equal(exact.scfIntegrals.type,'exact');assert.ok(exact.eri.length===exact.nao**4);assert.ok((exact.uhf??exact.rhf).converged);
}
for(const options of [{scfIntegrals:'bad'},{scfIntegrals:'df',jkfitBasis:'cc-pvdz-ri'}])assert.throws(()=>validateElectronicOptions(water,options));
for(const method of ['rhf','uhf']){
  const noDIIS=await computeElectronic(engine,water,{method,diis:false});assert.ok((noDIIS.rhf??noDIIS.uhf).converged);
  const limited=await computeElectronic(engine,water,{method,maxIterations:1});assert.equal((limited.rhf??limited.uhf).converged,false);
  let cancel=false;await assert.rejects(()=>computeElectronic(engine,water,{method,isCancelled:()=>cancel},()=>{cancel=true;}),{name:'AbortError'});
}
const {handle}=prepareIntegrals(engine,water,false);
try{assert.equal(engine._qc_data(handle,3),0);assert.equal(engine._qc_rhf_create(handle,10,100,1,1e-10,1e-8,1e-8),0);assert.match(engine.UTF8ToString(engine._qc_rhf_error()),/Exact SCF requires/);}finally{engine._qc_free(handle);}
assert.throws(()=>computeERI(engine,water,[-1,0,0,0]),/outside/);
assert.throws(()=>computeERI(engine,water,[0,0]),/four/);
for(const method of ['rhf','df-rmp2','df-ump2']){
  const scan=await runScan(engine,h2,{a:0,b:1,start:.6,end:.8,count:3},{method,guess:'core'});
  assert.equal(scan.spec.scfIntegrals,'df');assert.ok(scan.points.every(p=>p.converged&&p.result.scfIntegrals.type==='df'&&!p.result.eri));
  assert.match(scanCSV(scan),/scf_integrals,scf_fitting_basis/);assert.match(scanCSV(scan),/def2-universal-JKFIT/);assert.match(scanXYZ(scan),/scf_integrals=df jkfit_basis=def2-universal-JKFIT/);
  assert.match(JSON.stringify(pesFigureModel({...scan,source:h2})),/Density fitted \(DF-JK\).*def2-universal-JKFIT/);
}
console.log('Default DF, independent exact backend for all methods, no-DIIS, failed references, cancellation/recovery, validation, PES and figure metadata PASS');
