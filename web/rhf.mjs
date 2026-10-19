// SPDX-License-Identifier: AGPL-3.0-only
import { prepareIntegrals } from './engine.mjs';
export const RHF_THRESHOLDS = { energy:1e-10, density:1e-8, residual:1e-8 };

export async function computeRHF(engine, molecule, options = {}, onIteration = () => {}) {
  const charge = options.charge ?? 0, maxIterations = options.maxIterations ?? 100;
  if (!Number.isSafeInteger(charge)) throw new Error('Molecular charge must be a safe integer.');
  if (!Number.isInteger(maxIterations) || maxIterations < 1 || maxIterations > 500) throw new Error('Maximum iterations must be an integer from 1 to 500.');
  const electrons = molecule.atoms.reduce((sum,atom) => sum+atom.charge,0)-charge;
  if (!Number.isSafeInteger(electrons) || electrons < 2 || electrons % 2) throw new Error('RHF requires a positive even electron count (closed-shell singlet). Check the charge.');
  if (electrons > 80) throw new Error('Electron count exceeds the capacity of the 40 AO prototype.');
  const { handle, result } = prepareIntegrals(engine,molecule);
  let rhf = 0;
  try {
    if (!result.checks.passed) throw new Error('Integral normalization or symmetry checks failed. RHF was not started.');
    const t = RHF_THRESHOLDS;
    rhf = engine._qc_rhf_create(handle,electrons,maxIterations,options.diis === false ? 0 : 1,t.energy,t.density,t.residual);
    if (!rhf) throw new Error(engine.UTF8ToString(engine._qc_rhf_error()));
    const history = [];
    let status = 0;
    while (status === 0) {
      status = engine._qc_rhf_step(rhf);
      if (status < 0) throw new Error(engine.UTF8ToString(engine._qc_rhf_error()));
      const entry = { iteration:engine._qc_rhf_iteration(rhf), energy:engine._qc_rhf_scalar(rhf,0),
        deltaEnergy:engine._qc_rhf_scalar(rhf,1), densityRMS:engine._qc_rhf_scalar(rhf,2), residualRMS:engine._qc_rhf_scalar(rhf,3) };
      history.push(entry); onIteration(entry);
      // Yield between steps so progress is delivered and the worker remains responsive.
      await new Promise(resolve => setTimeout(resolve,0));
    }
    const n = result.nao;
    result.rhf = { method:'RHF', charge, electrons, spin:0, converged:status === 1,
      status:status === 1 ? 'converged' : 'not_converged', iterations:history.length,
      energy:engine._qc_rhf_scalar(rhf,0), electronicEnergy:engine._qc_rhf_scalar(rhf,0)-result.nuclearRepulsion,
      smallestOverlapEigenvalue:engine._qc_rhf_scalar(rhf,4), thresholds:t,
      diis:options.diis !== false, history };
    ['density','fock','coefficients','orbitalEnergies','occupations'].forEach((key,kind) => {
      const pointer = engine._qc_rhf_data(rhf,kind), size = kind < 3 ? n*n : n;
      result.rhf[key] = engine.HEAPF64.slice(pointer/8,pointer/8+size);
    });
    return result;
  } finally {
    if (rhf) engine._qc_rhf_free(rhf);
    engine._qc_free(handle);
  }
}
