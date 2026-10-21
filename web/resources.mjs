// SPDX-License-Identifier: AGPL-3.0-only
import {basisSummary,packBasis} from './molecule.mjs';
import {AUXILIARY_BASIS} from './auxiliary-basis.mjs';
import {JKFIT_BASIS} from './jkfit-basis.mjs';
import {METHOD_SPECS,isUnrestricted} from './methods.mjs';
import {MAX_AO,MAX_SHELLS,MAX_AUXILIARY,MIB,MAX_JOB_BYTES,DEFAULT_JOB_BYTES} from './limits.mjs';
export function deviceBudget(navigator=globalThis.navigator) {
  const gib=Number(navigator?.deviceMemory);
  // A coarse RAM hint, never a measurement of free memory. No active allocation probes.
  return Number.isFinite(gib)&&gib>0?Math.min(MAX_JOB_BYTES,Math.max(64*MIB,gib*64*MIB)):DEFAULT_JOB_BYTES;
}
export function jobBudget(value) {
  return Number.isFinite(value)&&value>0?Math.max(64*MIB,Math.min(MAX_JOB_BYTES,value)):DEFAULT_JOB_BYTES;
}
export const memoryText=bytes=>`${Math.ceil(bytes/MIB)} MiB`;
function auxiliarySize(molecule,catalog) {
  let nao=0,shells=0;
  for(const atom of molecule.atoms)for(const s of catalog.elements[atom.symbol].shells){nao+=2*s.l+1;++shells;}
  return {nao,shells};
}
export function estimateResources(molecule,options={}) {
  const basis=basisSummary(molecule),n=basis.nao,n2=n*n,integrals=options.mode==='integrals';
  const method=METHOD_SPECS[options.method??'rhf']??METHOD_SPECS.rhf,df=!integrals&&(options.scfIntegrals??'df')==='df';
  const jk=df?auxiliarySize(molecule,JKFIT_BASIS):{nao:0,shells:0};
  const ri=!integrals&&method.densityFitted?auxiliarySize(molecule,AUXILIARY_BASIS):{nao:0,shells:0};
  const steps=options.mode==='scan'&&Number.isInteger(options.scanCount)?Math.max(2,Math.min(101,options.scanCount)):1;
  const cycles=integrals?1:Number.isInteger(options.maxIterations)?Math.max(1,Math.min(500,options.maxIterations)):100;
  const electrons=molecule.atoms.reduce((a,b)=>a+b.charge,0)-(Number.isFinite(options.charge)?options.charge:0);
  const spin=isUnrestricted(options.method)?Math.max(0,(options.multiplicity??1)-1):0;
  const oa=Math.max(0,Math.min(n,(electrons+spin)/2)),ob=Math.max(0,Math.min(n,(electrons-spin)/2));
  const occ=isUnrestricted(options.method)?oa+ob:oa;
  // Upper bound uses untruncated metric rank. Include live copies and eigensolver workspace.
  const fitting=(m,multiplier=2)=>8*(multiplier*n2*m+6*m*m+48*n2);
  const denseBytes=8*n**4,referenceBytes=df?8*n2*jk.nao:2*denseBytes;
  let working=df?fitting(jk.nao):2*denseBytes+8*48*n2;
  if(ri.nao)working=Math.max(working,referenceBytes+fitting(ri.nao,3));
  if(!integrals&&options.method==='mp2')working=Math.max(working,referenceBytes+8*(3*n**3+48*n2));
  const resultBytes=8*20*n2+cycles*256+8192;
  // Worker and UI each retain streamed scan results. A scan runs sequentially.
  const scanBytes=steps>1?2*steps*resultBytes:resultBytes;
  const budgetBytes=jobBudget(options.resourceBudgetBytes),historyReserveBytes=Math.min(32*MIB,budgetBytes/8);
  const peakBytes=Math.ceil(24*MIB+historyReserveBytes+1.35*(working+scanBytes));
  const reasons=[],suggestions=[];
  if(n>MAX_AO)reasons.push(`This build supports at most ${MAX_AO} AOs.`);
  if(basis.nshell>MAX_SHELLS)reasons.push(`This basis exceeds ${MAX_SHELLS} orbital shells.`);
  if(packBasis(molecule,{enforceLimit:false}).exponents.length>4096)reasons.push('This basis exceeds 4096 orbital primitives.');
  for(const [label,size] of [['SCF',jk],['MP2',ri]])if(size.nao>MAX_AUXILIARY||size.shells>1024)reasons.push(`${label} fitting needs ${size.nao} auxiliary functions; this build supports ${MAX_AUXILIARY}.`);
  if(!df&&denseBytes>256*MIB)reasons.push('The exact ERI tensor exceeds the engine allocation bound.');
  if(df&&fitting(jk.nao)>384*MIB||ri.nao&&fitting(ri.nao)>384*MIB)reasons.push('The fitting workspace exceeds the engine allocation bound.');
  if(peakBytes>budgetBytes)reasons.push(`Estimated peak ${memoryText(peakBytes)} exceeds this device's ${memoryText(budgetBytes)} calculation budget.`);
  if(integrals)suggestions.push('Switch to a single-point DF-JK calculation if the full ERI tensor is unnecessary');
  else if(!df)suggestions.push('Use DF-JK SCF');
  if(!integrals&&options.method==='mp2')suggestions.push('use DF-RMP2 correlation');
  if(steps>2)suggestions.push('reduce scan points');
  suggestions.push('choose a smaller basis or molecule');
  // Relative effort, not seconds. Metric factorization and exchange remain expensive.
  const scfWork=df?jk.nao**3+n2*jk.nao*Math.max(1,oa+ob)*cycles:n**4*cycles;
  const correlation=ri.nao?ri.nao**3+n2*ri.nao*Math.max(1,occ)+Math.max(1,occ)**2*n2*ri.nao:!integrals&&options.method==='mp2'?Math.max(1,occ)*n**4:0;
  const workScore=(scfWork+correlation)*steps;
  const heavy=workScore>2e9||peakBytes>budgetBytes*.65, elevated=workScore>2e8||peakBytes>budgetBytes*.35;
  const level=reasons.length?'blocked':heavy?'heavy':elevated?'moderate':'light';
  const label={blocked:'Over budget',heavy:'Heavy calculation',moderate:'Moderate calculation',light:'Light calculation'}[level];
  const advice=suggestions.join('; ')+'.';
  return {nao:n,orbitalShells:basis.nshell,jkfitFunctions:jk.nao,correlationFunctions:ri.nao,steps,budgetBytes,peakBytes,workingBytes:working,resultBytes,historyReserveBytes,workScore,level,label,blocked:!!reasons.length,reasons,advice,
    summary:`${label} · est. ${memoryText(peakBytes)} peak${steps>1?` · ${steps} sequential points`:''}`,
    detail:`${reasons.join(' ')}${reasons.length?' ':''}${level==='heavy'?'May take a long time; Stop can restart the engine. ':''}${level==='blocked'||level==='heavy'?advice+' ':''}Allocation estimate, not free RAM or a runtime prediction. Other tabs and browser limits can still affect the calculation.`};
}
export function assertResources(molecule,options={}) {
  const estimate=estimateResources(molecule,options);
  if(estimate.blocked)throw new Error(estimate.reasons.join(' ')+' '+estimate.advice);
  return estimate;
}
// Count unique typed-array buffers in retained results, without serializing them.
export function retainedBytes(value,seen=new Set()) {
  if(!value||typeof value!=='object'||seen.has(value))return 0;
  seen.add(value);
  if(ArrayBuffer.isView(value)){if(seen.has(value.buffer))return 0;seen.add(value.buffer);return value.buffer.byteLength;}
  return Object.values(value).reduce((sum,item)=>sum+retainedBytes(item,seen),0);
}
