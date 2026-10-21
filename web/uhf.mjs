// SPDX-License-Identifier: AGPL-3.0-only
import {computeDFMP2} from './df-mp2.mjs';
import {prepareSCF} from './scf-integrals.mjs';
import {basisSummary} from './molecule.mjs';
import {RHF_THRESHOLDS,CalculationCancelled} from './rhf.mjs';
export function validateUHFOptions(molecule,options={}) {
  const charge=options.charge??0,maxIterations=options.maxIterations??100,multiplicity=options.multiplicity??1;
  if(!Number.isSafeInteger(charge)) throw new Error('Molecular charge must be a safe integer.');
  if(!Number.isInteger(maxIterations)||maxIterations<1||maxIterations>500) throw new Error('Maximum iterations must be an integer from 1 to 500.');
  const electrons=molecule.atoms.reduce((sum,a)=>sum+a.charge,0)-charge,spin=multiplicity-1;
  if(!Number.isSafeInteger(electrons)||electrons<1) throw new Error('UHF requires a positive electron count. Check the charge.');
  if(!Number.isInteger(multiplicity)||multiplicity<1||spin>electrons||(electrons+spin)%2) throw new Error('Multiplicity must be a positive integer consistent with the electron count: even electrons need odd multiplicity, odd electrons need even multiplicity.');
  const nalpha=(electrons+spin)/2,nbeta=(electrons-spin)/2;
  if(nalpha>basisSummary(molecule).nao) throw new Error('The spin populations do not fit this basis.');
  const guess=options.guess??'mix';
  if(!['core','mix'].includes(guess)) throw new Error('UHF guess must be core or mix.');
  return {charge,maxIterations,multiplicity,electrons,spin,nalpha,nbeta,guess};
}
export async function computeUHF(engine,molecule,options={},onIteration=()=>{},onProgress=()=>{}) {
  const check=()=>{if(options.isCancelled?.()) throw new CalculationCancelled();};check();
  const settings=validateUHFOptions(molecule,options),prepared=prepareSCF(engine,molecule,options),{handle,result}=prepared;let uhf=0,df=prepared.df;
  try {
    if(!result.checks.passed) throw new Error('Integral normalization or symmetry checks failed. UHF was not started.');
    const t=RHF_THRESHOLDS;
    uhf=engine._qc_uhf_create_df(handle,df,settings.nalpha,settings.nbeta,settings.maxIterations,options.diis===false?0:1,settings.guess==='mix'?1:0,t.energy,t.density,t.residual);
    if(!uhf) throw new Error(engine.UTF8ToString(engine._qc_uhf_error()));
    if(df){result.scfIntegrals.metricRank=engine._qc_uhf_df_rank(uhf);engine._qc_df_free(df);df=0;}
    const history=[];let status=0;
    while(status===0) {
      check();status=engine._qc_uhf_step(uhf);
      if(status<0) throw new Error(engine.UTF8ToString(engine._qc_uhf_error()));
      const entry={iteration:engine._qc_uhf_iteration(uhf),energy:engine._qc_uhf_scalar(uhf,0),deltaEnergy:engine._qc_uhf_scalar(uhf,1),densityRMS:engine._qc_uhf_scalar(uhf,2),residualRMS:engine._qc_uhf_scalar(uhf,3)};
      history.push(entry);onIteration(entry);await new Promise(resolve=>setTimeout(resolve,0));
    }
    check();const n=result.nao,value=i=>engine._qc_uhf_scalar(uhf,i);
    const r=result.uhf={method:'UHF',...settings,converged:status===1,status:status===1?'converged':'not_converged',iterations:history.length,energy:value(0),electronicEnergy:value(0)-result.nuclearRepulsion,smallestOverlapEigenvalue:value(4),spinSquare:value(5),targetSpinSquare:value(6),spinContamination:Math.max(0,value(5)-value(6)),mixedGuessApplied:!!value(7),thresholds:t,diis:options.diis!==false,history};
    ['densityAlpha','densityBeta','fockAlpha','fockBeta','coefficientsAlpha','coefficientsBeta','orbitalEnergiesAlpha','orbitalEnergiesBeta','occupationsAlpha','occupationsBeta'].forEach((key,kind)=>{
      const pointer=engine._qc_uhf_data(uhf,kind),size=kind<6?n*n:n;r[key]=engine.HEAPF64.slice(pointer/8,pointer/8+size);
    });
    r.density=r.densityAlpha.map((v,i)=>v+r.densityBeta[i]);r.spinDensity=r.densityAlpha.map((v,i)=>v-r.densityBeta[i]);
    if(options.dfmp2)result.mp2=r.converged?await computeDFMP2(engine,uhf,molecule,{...options,reference:'UHF'},onProgress):{method:'DF-UMP2',reference:'UHF',allElectron:true,densityFitted:true,auxiliaryBasis:'cc-pVDZ-RI',status:'not_computed',energy:null,reason:'UHF did not converge.'};
    return result;
  } finally {if(uhf) engine._qc_uhf_free(uhf);if(df)engine._qc_df_free(df);engine._qc_free(handle);}
}
