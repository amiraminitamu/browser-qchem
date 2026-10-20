// SPDX-License-Identifier: AGPL-3.0-only
// Geometry is always in bohr. Angular coordinates are signed degrees.
import {ANGSTROM_PER_BOHR} from './molecule.mjs';
export const COORDINATES=Object.freeze({
  bond:{label:'Bond length',symbol:'R',arity:2,unit:'Å',presets:[['Equilibrium bracket',null,null,.25],['Bond dissociation',.70,2.80]]},
  angle:{label:'Valence angle',symbol:'θ',arity:3,unit:'°',presets:[['Valence bend',null,null,25],['Linearization',90,180]]},
  dihedral:{label:'Dihedral angle',symbol:'φ',arity:4,unit:'°',presets:[['Full internal rotation',0,360],['Cis–trans isomerization',0,180]]}
});
export const sub=(a,b)=>a.map((x,i)=>x-b[i]);
export const add=(a,b)=>a.map((x,i)=>x+b[i]);
export const scale=(a,k)=>a.map(x=>x*k);
export const dot=(a,b)=>a.reduce((s,x,i)=>s+x*b[i],0);
export const cross=(a,b)=>[a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]];
export function normalize(a) {const n=Math.hypot(...a);if(n<1e-10)throw new Error('Coordinate is undefined for coincident or collinear atoms.');return scale(a,1/n);}
export function rotateVector(v,axis,angle) {
  return add(add(scale(v,Math.cos(angle)),scale(cross(axis,v),Math.sin(angle))),scale(axis,dot(axis,v)*(1-Math.cos(angle))));
}
export function coordinateValue(molecule,type,indices,unit='angstrom') {
  const [a,b,c,d]=indices.map(i=>molecule.atoms[i]?.xyz);
  if(!a||!b)throw new Error('Choose valid atoms for the coordinate.');
  if(type==='bond')return Math.hypot(...sub(b,a))*(unit==='bohr'?1:ANGSTROM_PER_BOHR);
  if(!c)throw new Error('Choose a third atom.');
  if(type==='angle'){const u=normalize(sub(a,b)),v=normalize(sub(c,b));return Math.atan2(Math.hypot(...cross(u,v)),dot(u,v))*180/Math.PI;}
  if(!d)throw new Error('Choose a fourth atom.');
  const axis=normalize(sub(c,b)),u=sub(sub(a,b),scale(axis,dot(sub(a,b),axis))),v=sub(sub(d,c),scale(axis,dot(sub(d,c),axis)));
  normalize(u);normalize(v);
  return Math.atan2(dot(axis,cross(u,v)),dot(u,v))*180/Math.PI;
}
export function displaceCoordinate(molecule,spec,value) {
  const type=spec.coordinate??'bond',indices=spec.atoms??[spec.a,spec.b];
  const initial=coordinateValue(molecule,type,indices,spec.unit);
  const xyz=indices.map(i=>molecule.atoms[i].xyz);
  let transform;
  if(type==='bond') {
    const shift=scale(normalize(sub(xyz[1],xyz[0])),(value-initial)/(spec.unit==='angstrom'?ANGSTROM_PER_BOHR:1));
    transform=position=>add(position,shift);
  } else {
    let axis;
    if(type==='dihedral')axis=normalize(sub(xyz[2],xyz[1]));
    else {
      axis=cross(sub(xyz[0],xyz[1]),sub(xyz[2],xyz[1]));
      // A linear bend has no unique plane; choose a deterministic perpendicular.
      if(Math.hypot(...axis)<1e-10){const arm=normalize(sub(xyz[0],xyz[1]));axis=cross(arm,Math.abs(arm[0])<.9?[1,0,0]:[0,1,0]);}
      axis=normalize(axis);
    }
    const angle=(value-initial)*Math.PI/180,origin=xyz[1];
    transform=position=>add(origin,rotateVector(sub(position,origin),axis,angle));
  }
  return {...molecule,atoms:molecule.atoms.map((atom,i)=>({...atom,xyz:spec.moving.includes(i)?transform(atom.xyz):[...atom.xyz]}))};
}
export function coordinateLabel(molecule,spec) {
  const type=spec.coordinate??'bond',meta=COORDINATES[type];
  const atoms=(spec.atoms??[spec.a,spec.b]).map(i=>`${molecule.atoms[i].symbol}${String(i+1).replace(/\d/g,d=>'₀₁₂₃₄₅₆₇₈₉'[d])}`).join('–');
  return `${meta.label} ${meta.symbol}(${atoms}) [${type==='bond'?(spec.unit==='bohr'?'bohr':'Å'):'°'}]`;
}
