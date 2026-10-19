// SPDX-License-Identifier: AGPL-3.0-only
import { WATER_XYZ } from './molecule.mjs';
import { serializable } from './engine.mjs';
const $ = id => document.getElementById(id);
let current = null, busy = true;
$('geometry').value = WATER_XYZ;
const worker = new Worker(new URL('./worker.mjs', import.meta.url), { type:'module' });
function setBusy(value) {
  busy = value; $('calculate').disabled = value; $('water').disabled = value;
  $('download').disabled = value || !current; $('geometry').disabled = value; $('unit').disabled = value;
  for (const id of ['mode','charge','maxiter','diis']) $(id).disabled = value;
}
function startHistory() {
  const table = document.createElement('table'), row = table.createTHead().insertRow();
  for (const title of ['Cycle','Total E / Eh','delta E / Eh','Density RMS','Residual RMS']) {
    const cell = document.createElement('th'); cell.textContent = title; row.append(cell);
  }
  table.createTBody(); $('scf-history').replaceChildren(table);
}
function appendIteration(entry) {
  const row = $('scf-history').querySelector('tbody').insertRow();
  for (const value of [entry.iteration,entry.energy.toFixed(12),entry.deltaEnergy.toExponential(3),entry.densityRMS.toExponential(3),entry.residualRMS.toExponential(3)]) row.insertCell().textContent = value;
  $('scf-history').scrollTop = $('scf-history').scrollHeight;
}
function showOrbitals(rhf) {
  const table = document.createElement('table'), header = table.insertRow();
  for (const title of ['Orbital','Energy / Eh','Occupation']) {
    const th = document.createElement('th'); th.textContent = title; header.append(th);
  }
  rhf.orbitalEnergies.forEach((energy,i) => {
    const row = table.insertRow();
    for (const value of [i+1,energy.toFixed(10),rhf.occupations[i]]) row.insertCell().textContent = value;
  });
  $('orbital-table').replaceChildren(table);
}
function run() {
  if (busy) return;
  if ($('mode').value === 'rhf' && (!$('charge').reportValidity() || !$('maxiter').reportValidity())) return;
  current = null; $('results').hidden = true;
  $('scf-results').hidden = $('mode').value !== 'rhf';
  startHistory(); $('orbital-table').textContent = 'Available after the SCF calculation.';
  $('status').textContent = 'Computing molecular integrals in the browser worker...';
  $('summary').textContent = 'Evaluating all one- and two-electron integrals.';
  setBusy(true); worker.postMessage({ geometry:$('geometry').value, unit:$('unit').value,
    mode:$('mode').value,options:{ charge:Number($('charge').value),maxIterations:Number($('maxiter').value),diis:$('diis').checked } });
}
function showMatrix() {
  if (!current) return;
  const { nao:n, aoLabels } = current, kind = $('matrix-kind').value;
  const values = kind === 'core' ? current.kinetic.map((v,i) => v+current.nuclear[i])
    : ['density','fock'].includes(kind) ? current.rhf?.[kind] : current[kind];
  if (!values) { $('matrix').textContent = 'Run RHF to inspect this matrix.'; return; }
  const table = document.createElement('table'), header = table.insertRow();
  const corner = document.createElement('th'); corner.textContent = 'AO'; header.append(corner);
  for (let j = 0; j < n; ++j) { const th = document.createElement('th'); th.textContent = j+1; header.append(th); }
  for (let i = 0; i < n; ++i) {
    const row = table.insertRow(), label = document.createElement('th');
    label.textContent = `${i+1}: ${aoLabels[i]}`; row.append(label);
    for (let j = 0; j < n; ++j) row.insertCell().textContent = values[i*n+j].toFixed(8);
  }
  $('matrix').replaceChildren(table);
}
function showERI() {
  if (!current) return;
  const inputs = [...$('eri-selectors').querySelectorAll('input')];
  if (inputs.some(input => !input.checkValidity() || !Number.isInteger(Number(input.value)))) {
    $('eri-value').textContent = 'Choose four valid integer AO indices.'; return;
  }
  const [i,j,k,l] = inputs.map(input => Number(input.value)-1), n = current.nao;
  $('eri-value').textContent = `(${i+1} ${j+1}|${k+1} ${l+1}) = ${current.eri[((i*n+j)*n+k)*n+l].toFixed(14)} Eh`;
}
function showResult(result) {
  current = result; const { nao:n, checks } = result;
  $('status').textContent = result.rhf
    ? result.rhf.converged ? 'CONVERGED: RHF/STO-3G.' : 'NOT CONVERGED: maximum iterations reached.'
    : checks.passed ? 'PASS: normalization and symmetry checks.' : 'FAIL: normalization or symmetry checks.';
  $('summary').textContent = [
    `Atoms                 ${result.molecule.atoms.length}`, `Basis                 STO-3G / real spherical`,
    `Shells / AOs          ${result.nshell} / ${n}`, `One-electron arrays   3 x ${n} x ${n}`,
    `ERI tensor            ${n} x ${n} x ${n} x ${n} (${n**4} values)`,
    `Nuclear repulsion     ${result.nuclearRepulsion.toFixed(12)} Eh`, `Elapsed               ${result.elapsedMs.toFixed(1)} ms`, '',
    `Max |Sii - 1|         ${checks.diagonalError.toExponential(2)}`, `Matrix symmetry       ${checks.matrixSymmetryError.toExponential(2)}`,
    `ERI symmetry          ${checks.eriSymmetryError.toExponential(2)}`,
  ].join('\n');
  $('scf-results').hidden = !result.rhf;
  if (result.rhf) {
    const r = result.rhf;
    $('summary').textContent = [r.converged ? `Total RHF energy      ${r.energy.toFixed(12)} Eh` : `Unconverged energy    ${r.energy.toFixed(12)} Eh`,
      `Electronic energy     ${r.electronicEnergy.toFixed(12)} Eh`,
      `Nuclear repulsion     ${result.nuclearRepulsion.toFixed(12)} Eh`,
      `Electrons / charge    ${r.electrons} / ${r.charge}`, `Iterations / DIIS     ${r.iterations} / ${r.diis ? 'on' : 'off'}`,
      `Basis / AOs           STO-3G / ${n}`,`Elapsed               ${result.elapsedMs.toFixed(1)} ms`,'',
      `Density RMS change    ${r.history.at(-1).densityRMS.toExponential(3)}`,
      `SCF residual RMS      ${r.history.at(-1).residualRMS.toExponential(3)}`,
      `Integral checks       ${checks.passed ? 'PASS' : 'FAIL'}`].join('\n');
    startHistory(); r.history.forEach(appendIteration); showOrbitals(r);
  }
  $('results').hidden = false;
  $('ao-order').textContent = result.aoLabels.map((label,i) => `${i+1}: ${label}`).join('\n');
  $('eri-selectors').replaceChildren();
  for (const name of ['i','j','k','l']) {
    const label = document.createElement('label'); label.textContent = name;
    const input = document.createElement('input');
    input.type = 'number'; input.min = '1'; input.max = n; input.step = '1'; input.value = '1'; input.required = true;
    input.setAttribute('aria-label',`ERI index ${name}`); input.addEventListener('input',showERI);
    label.append(input); $('eri-selectors').append(label);
  }
  showMatrix(); showERI();
}
worker.onmessage = ({ data }) => {
  if (data.type === 'ready') { setBusy(false); run(); return; }
  if (data.type === 'result') { showResult(data.result); setBusy(false); return; }
  if (data.type === 'iteration') {
    appendIteration(data.entry);
    $('status').textContent = `SCF iteration ${data.entry.iteration}: E = ${data.entry.energy.toFixed(10)} Eh`;
    $('summary').textContent = `SCF running...\nTotal energy     ${data.entry.energy.toFixed(12)} Eh\nDensity RMS      ${data.entry.densityRMS.toExponential(3)}\nResidual RMS     ${data.entry.residualRMS.toExponential(3)}`;
    return;
  }
  current = null; $('results').hidden = true;
  $('scf-results').hidden = true;
  $('status').textContent = data.type === 'fatal' ? 'Could not load the engine.' : 'Input or calculation error.';
  $('summary').textContent = data.message; setBusy(data.type === 'fatal');
};
worker.onerror = event => {
  current = null; $('results').hidden = true; $('status').textContent = 'The calculation worker failed.';
  $('scf-results').hidden = true;
  $('summary').textContent = `${event.message}\nCheck that the build completed and reload the page.`; setBusy(true);
};
$('calculate').addEventListener('click',run);
$('water').addEventListener('click',() => { $('geometry').value = WATER_XYZ; $('unit').value = 'angstrom'; $('charge').value = '0'; run(); });
$('matrix-kind').addEventListener('change',showMatrix);
for (const id of ['geometry','unit','mode','charge','maxiter','diis']) $(id).addEventListener('input',() => {
  current = null; $('results').hidden = true; $('download').disabled = true;
  $('scf-results').hidden = true;
  $('status').textContent = 'Input changed. Compute integrals to update the results.';
  $('summary').textContent = 'Results cleared because the input changed.';
});
$('download').addEventListener('click',() => {
  if (!current) return;
  const blob = new Blob([JSON.stringify(serializable(current),null,2)],{ type:'application/json' });
  const url = URL.createObjectURL(blob), link = document.createElement('a');
  link.href = url; link.download = current.rhf ? 'rhf-sto3g.json' : 'molecular-integrals-sto3g.json'; link.click(); setTimeout(() => URL.revokeObjectURL(url),1000);
});
setBusy(true);
