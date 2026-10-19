import assert from 'node:assert/strict';
import fs from 'node:fs';
import createQChem from '../build/wasm/dist/qc_core.mjs';
import { parseGeometry, packBasis, WATER_XYZ, ANGSTROM_PER_BOHR } from '../web/molecule.mjs';
import { computeIntegrals } from '../web/engine.mjs';

const engine = await createQChem();
const reference = JSON.parse(fs.readFileSync(new URL('./fixtures/water-sto3g.json',import.meta.url),'utf8'));
function compareWater() {
  const result = computeIntegrals(engine,parseGeometry(WATER_XYZ,'angstrom'));
  assert.equal(result.nao,7); assert.equal(result.eri.length,2401); assert.ok(result.checks.passed);
  for (const key of ['overlap','kinetic','nuclear','eri']) {
    assert.equal(result[key].length,reference[key].length);
    const error = result[key].reduce((max,value,i) => Math.max(max,Math.abs(value-reference[key][i])),0);
    assert.ok(error < 1e-10,`${key}: ${error}`);
    console.log(`Water ${key}: all ${result[key].length} entries, max error ${error.toExponential(3)} PASS`);
  }
  assert.ok(Math.abs(result.nuclearRepulsion-reference.nuclearRepulsion) < 1e-11);
  return result;
}
compareWater();
// Force memory growth, then recalculate to detect stale JavaScript heap views.
const before = engine.HEAPF64.buffer;
const large = engine._malloc(before.byteLength+1024);
assert.ok(large); assert.notEqual(engine.HEAPF64.buffer,before);
engine._free(large);
compareWater();
console.log('WASM memory growth and repeated calculation: PASS');
const bohr = parseGeometry('H 0 0 0\nH 0 0 1.4','bohr');
const angstrom = parseGeometry(`H 0 0 0\nH 0 0 ${1.4*ANGSTROM_PER_BOHR}`,'angstrom');
const h2 = computeIntegrals(engine,bohr), h2a = computeIntegrals(engine,angstrom);
assert.ok(Math.abs(h2.overlap[1]-engine._h2_overlap(1.4)) < 1e-12);
assert.deepEqual(h2.overlap,h2a.overlap);
for (const [text,unit] of [['Crap 0 0 0','bohr'],['H NaN 0 0','bohr'],['H 0 0','bohr'],
  ['2\ncomment\nH 0 0 0','bohr'],['H 0 0 0','meters']]) assert.throws(() => parseGeometry(text,unit));
assert.throws(() => computeIntegrals(engine,parseGeometry('H 0 0 0\nH 0 0 0','bohr')),/coincide/);
assert.throws(() => packBasis(parseGeometry(Array.from({length:9},(_,i) => `Ne ${i*3} 0 0`).join('\n'),'bohr')),/40 AOs/);
assert.equal(engine._qc_compute(0,0,0,0,0,0,0,0),0);
assert.match(engine.UTF8ToString(engine._qc_error()),/Invalid input/);
// Recover after rejected geometry, with no stale result handle.
compareWater();
console.log('Unit conversion, H2 regression, input rejection, and recovery: PASS');
