// SPDX-License-Identifier: AGPL-3.0-only
import {assertResources} from './resources.mjs';
import {scfBackend} from './scf-integrals.mjs';
import {METHOD_SPECS,isUnrestricted,isMP2} from './methods.mjs';
import {computeRHF,validateRHFOptions} from './rhf.mjs';
import {computeUHF,validateUHFOptions} from './uhf.mjs';
export function validateElectronicOptions(molecule,options={}) {
  assertResources(molecule,options);
  const method=options.method??'rhf',scfIntegrals=scfBackend(options);
  if(!Object.hasOwn(METHOD_SPECS,method))throw new Error('Choose RHF, UHF, MP2, DF-RMP2, or DF-UMP2.');
  if(METHOD_SPECS[method].densityFitted&&options.auxiliaryBasis&&options.auxiliaryBasis!=='cc-pvdz-ri')throw new Error('The supported correlation auxiliary basis is cc-pVDZ-RI.');
  try {return {...(isUnrestricted(method)?validateUHFOptions(molecule,options):validateRHFOptions(molecule,options)),method,scfIntegrals,...(METHOD_SPECS[method].densityFitted?{auxiliaryBasis:'cc-pvdz-ri'}:{})};}
  catch(error) {if(isMP2(method))throw new Error(method.toUpperCase()+' requires a valid '+METHOD_SPECS[method].reference+' reference: '+error.message);throw error;}
}
export async function computeElectronic(engine,molecule,options={},onIteration=()=>{},onProgress=()=>{}) {
  const {method}=validateElectronicOptions(molecule,options);
  const settings={...options,mp2:method==='mp2',dfmp2:!!METHOD_SPECS[method].densityFitted};
  const result=isUnrestricted(method)?await computeUHF(engine,molecule,settings,onIteration,onProgress):await computeRHF(engine,molecule,settings,onIteration,onProgress);
  result.method=method.toUpperCase();return result;
}
export const scfResult=result=>result?.uhf??result?.rhf;
export const totalEnergy=result=>result?.method==='MP2'||result?.mp2?result.mp2?.energy??null:scfResult(result)?.energy??null;
