// SPDX-License-Identifier: AGPL-3.0-only
export const SETTINGS_KEY='browser-qchem.settings.v2';
export const SETTINGS_FIELDS=['preset','geometry','unit','basis','method','mode','charge','multiplicity','guess','maxiter','diis','scan-coordinate','scan-a','scan-b','scan-c','scan-d','scan-moving','scan-start','scan-end','scan-count','scan-unit'];
const choices={preset:['water','h2','methane','ammonia','h2plus','heh','hehplus','custom'],unit:['angstrom','bohr'],basis:['sto-3g','6-31g','cc-pvdz'],method:['rhf','uhf','mp2'],mode:['rhf','scan','integrals'],guess:['core','mix'],'scan-unit':['angstrom','bohr','degree'],'scan-coordinate':['bond','angle','dihedral']};
const field=(form,id)=>form.querySelector('#'+id);
export function readSettings(form) {
  return Object.fromEntries(SETTINGS_FIELDS.map(id=>[id,id==='diis'?field(form,id).checked:field(form,id).value]));
}
export function saveSettings(form,storage,methodExplicit) {
  try {storage?.setItem(SETTINGS_KEY,JSON.stringify({schema:2,methodExplicit:!!methodExplicit,fields:readSettings(form)}));} catch {}
}
export function restoreSettings(form,storage) {
  try {
    const saved=JSON.parse(storage?.getItem(SETTINGS_KEY)??'null'),values=saved?.fields?{'scan-coordinate':'bond','scan-c':'','scan-d':'',...saved.fields}:null;
    if(saved?.schema!==2||!values||SETTINGS_FIELDS.some(id=>id==='diis'?typeof values[id]!=='boolean':typeof values[id]!=='string')) return false;
    if(Object.entries(choices).some(([id,allowed])=>!allowed.includes(values[id]))||values.geometry.length>1000000) return false;
    for(const id of SETTINGS_FIELDS) {if(id==='diis')field(form,id).checked=values[id];else field(form,id).value=values[id];}
    return !!saved.methodExplicit;
  } catch {return false;}
}
export function presetMethod(current,recommended,explicit) {return explicit?current:recommended??current;}
export function captureRequest(form,{id,protocol,buildId}) {
  const f=readSettings(form),method=f.method,mode=f.mode;
  if(!['rhf','uhf','mp2'].includes(method)) throw new Error('Choose RHF, UHF, or MP2.');
  if(!['rhf','scan','integrals'].includes(mode)) throw new Error('Choose a valid calculation type.');
  const options=Object.freeze({method,charge:Number(f.charge),multiplicity:Number(f.multiplicity),guess:f.guess,maxIterations:Number(f.maxiter),diis:f.diis});
  const moving=Object.freeze(f['scan-moving'].trim()?f['scan-moving'].trim().split(/[,\s]+/).map(v=>/^\d+$/.test(v)?Number(v)-1:NaN):[]);
  const scan=mode==='scan'?Object.freeze({coordinate:f['scan-coordinate'],atoms:Object.freeze(['a','b','c','d'].slice(0,({bond:2,angle:3,dihedral:4})[f['scan-coordinate']]??2).map(key=>f['scan-'+key]===''?NaN:Number(f['scan-'+key])-1)),a:Number(f['scan-a'])-1,b:Number(f['scan-b'])-1,moving,start:Number(f['scan-start']),end:Number(f['scan-end']),count:Number(f['scan-count']),unit:f['scan-unit']}):undefined;
  return Object.freeze({id,protocol,buildId,geometry:f.geometry,unit:f.unit,basis:f.basis,mode,options,scan});
}
export function assertResultMatches(request,result) {
  if(result.molecule?.basis!==request.basis) throw new Error('Engine returned a result for a different basis. Reload the page to restart the engine.');
  if(request.mode==='integrals') return;
  const expected=request.options.method.toUpperCase(),actual=result.method;
  if(actual!==expected) throw new Error('Engine returned '+(actual??'an unidentified method')+' for a '+expected+' request. Reload the page to restart the engine.');
  const scf=expected==='UHF'?result.uhf:result.rhf;
  if(!scf||scf.charge!==request.options.charge||(expected==='UHF'&&scf.multiplicity!==request.options.multiplicity)||(expected==='MP2'&&!result.mp2)) throw new Error('Engine result does not match the requested charge, spin or method. Reload the page to restart the engine.');
}
