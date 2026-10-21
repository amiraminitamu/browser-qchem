// SPDX-License-Identifier: AGPL-3.0-only
import {assertResources} from './resources.mjs?v=e933e699b7cf09943ce2';
import {prepareIntegrals,prepareDFIntegrals} from './engine.mjs?v=e933e699b7cf09943ce2';
import {packAuxiliary} from './auxiliary-basis.mjs?v=e933e699b7cf09943ce2';
import {JKFIT_BASIS} from './jkfit-basis.mjs?v=e933e699b7cf09943ce2';
export function scfBackend(options={}) {
  const type=options.scfIntegrals??'df';
  if(!['df','exact'].includes(type))throw new Error('SCF integrals must be density fitted or exact.');
  if(options.jkfitBasis&&options.jkfitBasis!==JKFIT_BASIS.id)throw new Error('The supported SCF fitting basis is def2-universal-JKFIT.');
  return type;
}
export function prepareSCF(engine,molecule,options={}) {
  const resources=assertResources(molecule,options);
  const type=scfBackend(options),prepared=prepareIntegrals(engine,molecule,type==='exact',options);let df=0;
  try {
    if(type==='df')df=prepareDFIntegrals(engine,molecule,packAuxiliary(molecule,JKFIT_BASIS));
    prepared.result.resources=resources;
    prepared.result.scfIntegrals={type,densityFitted:type==='df',...(df?{auxiliaryBasis:JKFIT_BASIS.label,auxiliaryFunctions:engine._qc_df_naux(df),metricEigenvalueCutoff:1e-10}:{}),storage:type==='df'?'three-index AO factors; no dense four-center ERI tensor':'dense four-center AO ERIs'};
    return {...prepared,df};
  }catch(error){if(df)engine._qc_df_free(df);engine._qc_free(prepared.handle);throw error;}
}
