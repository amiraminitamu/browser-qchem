// SPDX-License-Identifier: AGPL-3.0-only
import {AUXILIARY_BASIS,packAuxiliary} from './auxiliary-basis.mjs?v=e933e699b7cf09943ce2';
import {prepareDFIntegrals} from './engine.mjs?v=e933e699b7cf09943ce2';
export async function computeDFMP2(engine,reference,molecule,options={},onProgress=()=>{}) {
  const unrestricted=options.reference==='UHF',method=unrestricted?'DF-UMP2':'DF-RMP2';
  const check=()=>{if(options.isCancelled?.()){const error=new Error('Calculation cancelled.');error.name='AbortError';throw error;}};
  let df=0,handle=0;
  try {
    check();onProgress({phase:'fitting',method,completed:0,occupied:0});await new Promise(resolve=>setTimeout(resolve,0));check();
    df=prepareDFIntegrals(engine,molecule,packAuxiliary(molecule));
    check();onProgress({phase:'metric',method,completed:0,occupied:0});await new Promise(resolve=>setTimeout(resolve,0));check();
    handle=unrestricted?engine._qc_dfmp2_uhf_create(reference,df):engine._qc_dfmp2_rhf_create(reference,df);
    if(!handle)throw new Error(engine.UTF8ToString(engine._qc_dfmp2_error()));
    const auxiliaryFunctions=engine._qc_df_naux(df),metricRank=engine._qc_dfmp2_rank(handle),steps=engine._qc_dfmp2_steps(handle);
    engine._qc_df_free(df);df=0; // The post-SCF handle now owns its MO factors.
    let status=engine._qc_dfmp2_status(handle);
    while(status===0){
      check();status=engine._qc_dfmp2_step(handle);
      if(status<0)throw new Error(engine.UTF8ToString(engine._qc_dfmp2_error()));
      onProgress({phase:'correlation',method,completed:engine._qc_dfmp2_iteration(handle),occupied:steps,correlationEnergy:engine._qc_dfmp2_scalar(handle,1)});
      await new Promise(resolve=>setTimeout(resolve,0));
    }
    check();const value=kind=>{const number=engine._qc_dfmp2_scalar(handle,kind);return Number.isFinite(number)?number:null;};
    return {method,reference:unrestricted?'UHF':'RHF',allElectron:true,densityFitted:true,status:'complete',
      auxiliaryBasis:AUXILIARY_BASIS.label,auxiliaryFunctions,metricRank,metricEigenvalueCutoff:value(9),smallestMetricEigenvalue:value(10),
      energy:value(0),correlationEnergy:value(1),oppositeSpinEnergy:value(2),sameSpinEnergy:value(3),sameSpinAlphaEnergy:value(7),sameSpinBetaEnergy:value(8),
      minimumAbsoluteDenominator:value(4),referenceGap:value(5),referenceEnergy:value(6),occupied:steps,
      occupiedAlpha:engine._qc_dfmp2_nocc(handle,0),occupiedBeta:engine._qc_dfmp2_nocc(handle,1),virtualAlpha:engine._qc_dfmp2_nvirt(handle,0),virtualBeta:engine._qc_dfmp2_nvirt(handle,1)};
  }finally{if(handle)engine._qc_dfmp2_free(handle);if(df)engine._qc_df_free(df);}
}
