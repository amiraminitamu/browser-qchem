import assert from 'node:assert/strict';
import createQChem from '../build/wasm/dist/qc_core.mjs';

const engine = await createQChem();
const a = [3.42525091, 0.62391373, 0.16885540];
const c = [0.15432897, 0.53532814, 0.44463454];
function analytical(r) {
  let sum = 0, norm = 0;
  for (let i = 0; i < 3; ++i) {
    for (let j = 0; j < 3; ++j) {
      const w = c[i] * c[j] * (2 * Math.sqrt(a[i] * a[j]) / (a[i] + a[j])) ** 1.5;
      norm += w;
      sum += w * Math.exp(-a[i] * a[j] / (a[i] + a[j]) * r * r);
    }
  }
  return sum / norm;
}
for (const r of [0, 0.1, 0.7, 1.4, 2, 4, 10, 20]) {
  const value = engine._h2_overlap(r);
  assert.ok(Number.isFinite(value) && Math.abs(value - analytical(r)) < 1e-12,
    `overlap mismatch at ${r} bohr: ${value}`);
  console.log(`${r.toFixed(1)} bohr: ${value.toFixed(15)} PASS`);
}
for (const invalid of [-1, NaN, Infinity, -Infinity]) {
  assert.ok(Number.isNaN(engine._h2_overlap(invalid)));
}
console.log('Invalid inputs return NaN: PASS');
