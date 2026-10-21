// SPDX-License-Identifier: AGPL-3.0-only
import {computeDFMP2} from './df-mp2.mjs';
import { computeMP2 } from './mp2.mjs';
import {prepareSCF} from './scf-integrals.mjs';
import { basisSummary } from './molecule.mjs';
export const RHF_THRESHOLDS = { energy:1e-10, density:1e-8, residual:1e-8 };
export class CalculationCancelled extends Error {
  constructor() { super('Calculation cancelled.'); this.name='AbortError'; }
}

export function validateRHFOptions(molecule,options={}) {
  const charge = options.charge ?? 0, maxIterations = options.maxIterations ?? 100;
  if (!Number.isSafeInteger(charge)) throw new Error('Molecular charge must be a safe integer.');
  if (!Number.isInteger(maxIterations) || maxIterations < 1 || maxIterations > 500) throw new Error('Maximum iterations must be an integer from 1 to 500.');
  const electrons = molecule.atoms.reduce((sum,atom) => sum+atom.charge,0)-charge;
  if (!Number.isSafeInteger(electrons) || electrons < 2 || electrons % 2) throw new Error('RHF requires a positive even electron count (closed-shell singlet). Check the charge.');
  if (electrons > 2*basisSummary(molecule).nao) throw new Error('The electron count does not fit this basis. Check the charge.');
  return {charge,maxIterations,electrons};
}
export async function computeRHF(engine, molecule, options = {}, onIteration = () => {}, onProgress = () => {}) {
  const checkCancelled = () => { if (options.isCancelled?.()) throw new CalculationCancelled(); };
  checkCancelled();
  const {charge,maxIterations,electrons}=validateRHFOptions(molecule,options);
  const prepared=prepareSCF(engine,molecule,options),{handle,result}=prepared;
  let df=prepared.df;
  let rhf = 0;
  try {
    if (!result.checks.passed) throw new Error('Integral normalization or symmetry checks failed. RHF was not started.');
    const t = RHF_THRESHOLDS;
    rhf = engine._qc_rhf_create_df(handle,df,electrons,maxIterations,options.diis === false ? 0 : 1,t.energy,t.density,t.residual);
    if (!rhf) throw new Error(engine.UTF8ToString(engine._qc_rhf_error()));
    if(df){result.scfIntegrals.metricRank=engine._qc_rhf_df_rank(rhf);engine._qc_df_free(df);df=0;}
    const history = [];
    let status = 0;
    while (status === 0) {
      checkCancelled();
      status = engine._qc_rhf_step(rhf);
      if (status < 0) throw new Error(engine.UTF8ToString(engine._qc_rhf_error()));
      const entry = { iteration:engine._qc_rhf_iteration(rhf), energy:engine._qc_rhf_scalar(rhf,0),
        deltaEnergy:engine._qc_rhf_scalar(rhf,1), densityRMS:engine._qc_rhf_scalar(rhf,2), residualRMS:engine._qc_rhf_scalar(rhf,3) };
      history.push(entry); onIteration(entry);
      // Yield between steps so progress is delivered and the worker remains responsive.
      await new Promise(resolve => setTimeout(resolve,0));
    }
    checkCancelled();
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
    if(options.dfmp2)result.mp2=result.rhf.converged?await computeDFMP2(engine,rhf,molecule,{...options,reference:'RHF'},onProgress):{method:'DF-RMP2',reference:'RHF',allElectron:true,densityFitted:true,auxiliaryBasis:'cc-pVDZ-RI',status:'not_computed',energy:null,reason:'RHF did not converge.'};
    else if(options.mp2) result.mp2=result.rhf.converged?await computeMP2(engine,rhf,options,onProgress):{method:'MP2',reference:'RHF',allElectron:true,densityFitted:false,status:'not_computed',energy:null,reason:'RHF did not converge.'};
    return result;
  } finally {
    if (rhf) engine._qc_rhf_free(rhf);
    if(df)engine._qc_df_free(df);engine._qc_free(handle);
  }
}

