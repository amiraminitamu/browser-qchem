// SPDX-License-Identifier: AGPL-3.0-only
import {computeRHF,validateRHFOptions} from './rhf.mjs';
import {computeUHF,validateUHFOptions} from './uhf.mjs';
export function validateElectronicOptions(molecule,options={}) {
  const method=options.method??'rhf';
  if(!['rhf','uhf','mp2'].includes(method)) throw new Error('Choose RHF, UHF, or MP2.');
  try {return {...(method==='uhf'?validateUHFOptions(molecule,options):validateRHFOptions(molecule,options)),method};}
  catch(error) {if(method==='mp2')throw new Error('MP2 requires a valid closed-shell RHF reference: '+error.message);throw error;}
}
export async function computeElectronic(engine,molecule,options={},onIteration=()=>{},onProgress=()=>{}) {
  const {method}=validateElectronicOptions(molecule,options);
  const result=method==='uhf'?await computeUHF(engine,molecule,options,onIteration):await computeRHF(engine,molecule,{...options,mp2:method==='mp2'},onIteration,onProgress);
  result.method=method.toUpperCase();return result;
}
export const scfResult=result=>result?.uhf??result?.rhf;
export const totalEnergy=result=>result?.method==='MP2'||result?.mp2?result.mp2?.energy??null:scfResult(result)?.energy??null;
