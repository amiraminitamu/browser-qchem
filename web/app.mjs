import createQChem from './qc_core.mjs';

const status = document.querySelector('#status');
const result = document.querySelector('#result');
const button = document.querySelector('#calculate');
const distance = document.querySelector('#distance');

// Independent closed-form reference for normalized contracted s Gaussians.
function referenceOverlap(r) {
  const a = [3.42525091, 0.62391373, 0.16885540];
  const c = [0.15432897, 0.53532814, 0.44463454];
  let sum = 0, norm = 0;
  for (let i = 0; i < 3; ++i) {
    for (let j = 0; j < 3; ++j) {
      const weight = c[i] * c[j] * (2 * Math.sqrt(a[i] * a[j]) / (a[i] + a[j])) ** 1.5;
      norm += weight;
      sum += weight * Math.exp(-a[i] * a[j] / (a[i] + a[j]) * r * r);
    }
  }
  return sum / norm;
}

try {
  const engine = await createQChem();
  function calculate() {
    if (!distance.reportValidity()) return;
    const r = Number(distance.value);
    if (!Number.isFinite(r) || r < 0 || r > 20) {
      status.textContent = 'Enter a finite separation between 0 and 20 bohr.';
      return;
    }
    const overlap = engine._h2_overlap(r);
    const reference = referenceOverlap(r);
    const error = Math.abs(overlap - reference);
    const passed = Number.isFinite(overlap) && error < 1e-12;
    status.textContent = passed ? 'PASS: libcint agrees with the analytical reference.'
      : 'FAIL: integral differs from the analytical reference.';
    result.textContent = [
      `Separation       ${r.toFixed(6)} bohr`,
      `libcint / WASM   ${overlap.toFixed(15)}`,
      `Reference        ${reference.toFixed(15)}`,
      `Absolute error   ${error.toExponential(3)}`,
      '', 'S = [[1, S12], [S12, 1]]',
    ].join('\n');
  }
  button.disabled = false;
  document.querySelector('#integral-form').addEventListener('submit', event => {
    event.preventDefault();
    calculate();
  });
  calculate();
} catch (error) {
  status.textContent = 'Could not load the integral engine.';
  result.textContent = `${error.message}\n\nRun bash scripts/build_wasm.sh first, then serve build/wasm/dist over HTTP.`;
}
