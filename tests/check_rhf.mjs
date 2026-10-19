import assert from 'node:assert/strict';
import fs from 'node:fs';
import createQChem from '../build/wasm/dist/qc_core.mjs';
import {parseGeometry,WATER_XYZ} from '../web/molecule.mjs';
import {computeRHF} from '../web/rhf.mjs';
const engine=await createQChem();
const reference=JSON.parse(fs.readFileSync(new URL('./fixtures/water-rhf-sto3g.json',import.meta.url),'utf8'));
const water=parseGeometry(WATER_XYZ);
let progress=0;
const result=await computeRHF(engine,water,{},() => progress++);
const r=result.rhf;
assert.ok(r.converged); assert.equal(r.iterations,progress); assert.equal(r.electrons,10);
assert.ok(Math.abs(r.energy-reference.energy)<1e-9);
for (const key of ['density','fock','orbitalEnergies']) {
  const error=r[key].reduce((max,value,i) => Math.max(max,Math.abs(value-reference[key][i])),0);
  assert.ok(error<1e-7,`${key}: ${error}`); console.log(`Water ${key} max error ${error.toExponential(3)} PASS`);
}
const last=r.history.at(-1);
assert.ok(Math.abs(last.deltaEnergy)<r.thresholds.energy);
assert.ok(last.densityRMS<r.thresholds.density); assert.ok(last.residualRMS<r.thresholds.residual);
console.log(`Water total energy ${r.energy.toFixed(12)} Eh, ${r.iterations} iterations PASS`);
const limited=await computeRHF(engine,water,{maxIterations:1});
assert.equal(limited.rhf.converged,false); assert.equal(limited.rhf.status,'not_converged');
assert.equal(limited.rhf.iterations,1);
assert.ok((await computeRHF(engine,water,{diis:false})).rhf.converged);
for (const options of [{charge:1},{charge:.5},{charge:-1000},{maxIterations:0},{maxIterations:501}]) {
  await assert.rejects(() => computeRHF(engine,water,options));
}
await assert.rejects(() => computeRHF(engine,parseGeometry('H 0 0 0\nH 0 0 .00001','bohr')),/linearly dependent/);
const ion=await computeRHF(engine,parseGeometry('He 0 0 0\nH 0 0 1.5','bohr'),{charge:1});
assert.ok(ion.rhf.converged); assert.equal(ion.rhf.electrons,2);
assert.equal(engine._qc_rhf_step(0),-1);
assert.equal(engine._qc_rhf_create(0,2,100,1,1e-10,1e-8,1e-8),0);
assert.ok((await computeRHF(engine,water)).rhf.converged);
console.log('DIIS off, iteration limit, charge, invalid inputs, overlap rejection, recovery: PASS');
