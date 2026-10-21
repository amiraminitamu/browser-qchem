// SPDX-License-Identifier: AGPL-3.0-only
export async function computeMP2(engine,reference,options={},onProgress=()=>{}) {
  const check=()=>{if(options.isCancelled?.()) {const e=new Error('Calculation cancelled.');e.name='AbortError';throw e;}};
  check();
  const handle=engine._qc_mp2_create(reference);
  if(!handle) throw new Error(engine.UTF8ToString(engine._qc_mp2_error()));
  try {
    const occupied=engine._qc_mp2_nocc(handle),virtual=engine._qc_mp2_nvirt(handle);
    let status=engine._qc_mp2_status(handle);
    while(status===0) {
      check();status=engine._qc_mp2_step(handle);
      if(status<0) throw new Error(engine.UTF8ToString(engine._qc_mp2_error()));
      onProgress({completed:engine._qc_mp2_iteration(handle),occupied,correlationEnergy:engine._qc_mp2_scalar(handle,1)});
      await new Promise(resolve=>setTimeout(resolve,0));
    }
    check();
    const value=i=>{const v=engine._qc_mp2_scalar(handle,i);return Number.isFinite(v)?v:null;};
    return {method:'MP2',reference:'RHF',allElectron:true,densityFitted:false,status:'complete',
      energy:value(0),correlationEnergy:value(1),oppositeSpinEnergy:value(2),sameSpinEnergy:value(3),
      minimumAbsoluteDenominator:value(4),referenceGap:value(5),referenceEnergy:value(6),occupied,virtual};
  } finally {engine._qc_mp2_free(handle);}
}
