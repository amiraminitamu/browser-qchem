// SPDX-License-Identifier: AGPL-3.0-only
import { ANGSTROM_PER_BOHR, toXYZ, basisSummary } from './molecule.mjs';
import {COORDINATES,coordinateValue,displaceCoordinate} from './coordinates.mjs';
import {computeElectronic,validateElectronicOptions,scfResult,totalEnergy} from './electronic.mjs';
export const HARTREE_TO_KCAL = 627.5094740631;
const RADII = { H:.31,He:.28,Li:1.28,Be:.96,B:.84,C:.76,N:.71,O:.66,F:.57,Ne:.58 };
export function distance(molecule,a,b) {
  return Math.hypot(...molecule.atoms[a].xyz.map((x,d) => x-molecule.atoms[b].xyz[d]));
}
export function inferBonds(molecule) {
  const bonds=[];
  for (let a=0;a<molecule.atoms.length;a++) for (let b=0;b<a;b++) {
    const length=distance(molecule,a,b)*ANGSTROM_PER_BOHR;
    if (length>.1 && length<1.25*(RADII[molecule.atoms[a].symbol]+RADII[molecule.atoms[b].symbol])) bonds.push([a,b]);
  }
  return bonds;
}
export function connectedFragment(molecule,a,b) {
  if (a===b || !molecule.atoms[a] || !molecule.atoms[b]) throw new Error('Select two different atoms.');
  const edges=inferBonds(molecule).filter(([i,j]) => !((i===a&&j===b)||(i===b&&j===a)));
  const seen=new Set([b]), queue=[b];
  while (queue.length) {
    const next=queue.shift();
    for (const [i,j] of edges) {
      const other=i===next ? j : j===next ? i : -1;
      if (other>=0 && !seen.has(other)) { seen.add(other); queue.push(other); }
    }
  }
  if (seen.has(a)) throw new Error('The selected bond is in a ring or the fragments remain connected. Specify moving atoms explicitly.');
  return [...seen].sort((i,j)=>i-j);
}
export function scanSpec(molecule,config) {
  const coordinate=config.coordinate??'bond',meta=COORDINATES[coordinate];
  if(!meta)throw new Error('Choose a bond, angle, or dihedral coordinate.');
  const atoms=(config.atoms??[config.a,config.b]).map(Number),n=molecule.atoms.length;
  if(atoms.length!==meta.arity||new Set(atoms).size!==meta.arity||!atoms.every(i=>Number.isInteger(i)&&i>=0&&i<n))throw new Error(`Select ${meta.arity} different valid atoms in order.`);
  const [a,b]=atoms,terminal=atoms.at(-1),unit=config.unit??(coordinate==='bond'?'angstrom':'degree');
  if(!(coordinate==='bond'?['angstrom','bohr']:['degree']).includes(unit))throw new Error('Invalid coordinate unit.');
  const start=Number(config.start),end=Number(config.end),count=Number(config.count);
  if(![start,end].every(Number.isFinite)||start===end)throw new Error('Scan endpoints must be finite and different.');
  if(coordinate==='bond'&&![start,end].every(x=>x>0&&x<1e4))throw new Error('Scan endpoints must be different positive distances below 10000.');
  if(coordinate==='angle'&&![start,end].every(x=>x>0&&x<=180))throw new Error('Valence angles must be greater than 0° and at most 180°.');
  if(coordinate==='dihedral'&&![start,end].every(x=>Math.abs(x)<=3600))throw new Error('Torsion endpoints must lie within ±3600°.');
  if(!Number.isInteger(count)||count<2||count>101)throw new Error('Use between 2 and 101 scan points.');
  // Preserve the original bond API. The UI supplies connected fragments explicitly.
  const moving=[...new Set(config.moving??[terminal])],fixed=coordinate==='dihedral'?atoms.slice(0,2):atoms.slice(0,-1);
  if(!moving.length||!moving.every(i=>Number.isInteger(i)&&i>=0&&i<n)||!moving.includes(terminal)||fixed.some(i=>moving.includes(i)))throw new Error('Moving atoms must include the last selected atom and exclude fixed coordinate atoms.');
  coordinateValue(molecule,coordinate,atoms,unit);
  if(coordinate==='bond'&&distance(molecule,a,b)<1e-8)throw new Error('The selected atoms coincide.');
  if(!basisSummary(molecule).fits)throw new Error('This scan exceeds the 40 AO limit.');
  return {a,b,atoms,coordinate,start,end,count,unit,moving:moving.sort((i,j)=>i-j),basis:molecule.basis,kind:`unrelaxed ${coordinate} scan`};
}
export function scanGeometry(molecule,spec,value) {return displaceCoordinate(molecule,spec,value);}
export async function runScan(engine,molecule,config,options={},onPoint=()=>{},onIteration=()=>{},onProgress=()=>{}) {
  const settings=validateElectronicOptions(molecule,options);
  const spec={...scanSpec(molecule,config),method:settings.method,scfIntegrals:settings.scfIntegrals,jkfitBasis:settings.scfIntegrals==='df'?'def2-universal-JKFIT':null,charge:settings.charge,electrons:settings.electrons,multiplicity:settings.multiplicity??1,guess:settings.guess??'core',...(settings.auxiliaryBasis?{auxiliaryBasis:settings.auxiliaryBasis}:{})}, points=[];
  for (let index=0;index<spec.count;index++) {
    if (options.isCancelled?.()) break;
    const coordinate=spec.start+(spec.end-spec.start)*index/(spec.count-1);
    const geometry=scanGeometry(molecule,spec,coordinate);
    let point;
    try {
      const started=performance.now();
      // Independent guesses at every point, with the selected UHF mixing policy.
      const result=await computeElectronic(engine,geometry,options,entry=>onIteration({ index,entry }),progress=>onProgress({index,progress}));
      result.elapsedMs=performance.now()-started;
      const { eri,...compact }=result; // Do not retain N^4 arrays for an entire scan.
      point={ index,coordinate,...(spec.coordinate==='bond'?{coordinateBohr:coordinate/(spec.unit==='angstrom'?ANGSTROM_PER_BOHR:1)}:{coordinateDegrees:coordinate}),
        molecule:geometry,energy:totalEnergy(result),status:scfResult(result).status,converged:scfResult(result).converged,
        iterations:scfResult(result).iterations,result:compact };
    } catch (error) {
      if (error.name==='AbortError') break;
      point={ index,coordinate,...(spec.coordinate==='bond'?{coordinateBohr:coordinate/(spec.unit==='angstrom'?ANGSTROM_PER_BOHR:1)}:{coordinateDegrees:coordinate}),molecule:geometry,energy:null,status:'error',converged:false,iterations:0,error:error.message };
    }
    points.push(point); onPoint(point,spec);
    await new Promise(resolve=>setTimeout(resolve,0));
  }
  return { spec,points,cancelled:points.length<spec.count,requested:spec.count,completed:points.length };
}
export function scanCSV(scan) {
  const escape=value=>`"${String(value??'').replaceAll('"','""')}"`;
  return ['point,coordinate,coordinate_unit,energy_hartree,converged,status,iterations,basis,error,method,multiplicity,spin_square,correlation_energy_hartree,reference_energy_hartree,charge,electrons,coordinate_type,atom_indices,auxiliary_basis,auxiliary_functions,metric_rank,density_fitted,scf_integrals,scf_fitting_basis,scf_auxiliary_functions,scf_metric_rank',
    ...scan.points.map(p=>[p.index+1,p.coordinate,scan.spec.unit,p.energy??'',p.converged,p.status,p.iterations,scan.spec.basis,p.error??'',p.result?.method??(scan.spec.method??'rhf').toUpperCase(),p.result?.uhf?.multiplicity??scan.spec.multiplicity??1,p.result?.uhf?.spinSquare??'',p.result?.mp2?.correlationEnergy??'',p.result?.mp2?.referenceEnergy??'',p.result?.uhf?.charge??p.result?.rhf?.charge??scan.spec.charge??'',p.result?.uhf?.electrons??p.result?.rhf?.electrons??scan.spec.electrons??'',scan.spec.coordinate??'bond',(scan.spec.atoms??[scan.spec.a,scan.spec.b]).map(i=>i+1).join(' '),p.result?.mp2?.auxiliaryBasis??scan.spec.auxiliaryBasis??'',p.result?.mp2?.auxiliaryFunctions??'',p.result?.mp2?.metricRank??'',p.result?.mp2?.densityFitted??!!scan.spec.auxiliaryBasis,p.result?.scfIntegrals?.type??scan.spec.scfIntegrals,p.result?.scfIntegrals?.auxiliaryBasis??scan.spec.jkfitBasis??'',p.result?.scfIntegrals?.auxiliaryFunctions??'',p.result?.scfIntegrals?.metricRank??''].map(escape).join(','))].join('\n')+'\n';
}
export function scanXYZ(scan) {
  return scan.points.map(p=>toXYZ(p.molecule,'angstrom',`point=${p.index+1} coordinate=${scan.spec.coordinate??'bond'}:${p.coordinate} ${scan.spec.unit} energy=${p.energy??'NA'} Eh method=${p.result?.method??(scan.spec.method??'rhf').toUpperCase()} multiplicity=${p.result?.uhf?.multiplicity??scan.spec.multiplicity??1} charge=${p.result?.uhf?.charge??p.result?.rhf?.charge??scan.spec.charge??'NA'} electrons=${p.result?.uhf?.electrons??p.result?.rhf?.electrons??scan.spec.electrons??'NA'} scf_integrals=${scan.spec.scfIntegrals}${scan.spec.jkfitBasis?' jkfit_basis='+scan.spec.jkfitBasis:''}${scan.spec.auxiliaryBasis?' auxiliary_basis='+scan.spec.auxiliaryBasis:''} status=${p.status}`)).join('');
}
