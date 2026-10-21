// SPDX-License-Identifier: AGPL-3.0-only
import {ANGSTROM_PER_BOHR} from './molecule.mjs?v=e933e699b7cf09943ce2';
import {COORDINATES,coordinateValue} from './coordinates.mjs?v=e933e699b7cf09943ce2';
import {connectedFragment} from './scan.mjs?v=e933e699b7cf09943ce2';
const $=id=>document.getElementById(id);
export class PESSetup {
  constructor({getMolecule,isBusy,changed}) {
    this.getMolecule=getMolecule;this.isBusy=isBusy;this.changed=changed;
    for(const [order,id] of ['a','b','c','d'].entries()) {
      const token=document.createElement('label');token.className='atom-token';token.dataset.order=order;
      const number=document.createElement('span');number.textContent=`${order+1}:`;
      const symbol=document.createElement('span');symbol.className='atom-symbol';
      const input=document.createElement('input');Object.assign(input,{type:'number',id:`scan-${id}`,min:'1',max:'32',step:'1',value:order<2?String(order+1):''});
      input.defaultValue=input.value;
      input.setAttribute('aria-label',`Atom ${order+1} index`);input.addEventListener('change',()=>{if(this.complete())this.autoMoving();});
      token.append(number,symbol,input);$('atom-tokens').append(token);
    }
    document.querySelectorAll('[data-coordinate]').forEach(button=>button.addEventListener('click',()=>this.setKind(button.dataset.coordinate)));
    $('scan-clear').addEventListener('click',()=>this.clear());
    $('scan-reverse').addEventListener('click',()=>{if(this.isBusy()||!this.complete())return;this.setIndices(this.indices().reverse());this.autoMoving();this.changed();});
    $('scan-reset').addEventListener('click',()=>this.resetRange());
    $('scan-preset').addEventListener('change',()=>{const index=$('scan-preset').value;if(index==='')return;this.preset(Number(index));$('scan-preset').value='';});
    $('auto-fragment').addEventListener('click',()=>{if(this.isBusy())return;this.autoMoving();this.changed();});
    document.addEventListener('keydown',event=>{
      if(event.key!=='Escape'||$('mode').value!=='scan'||this.isBusy()||$('workspace-guide').open)return;
      if(event.target.closest('input,textarea,select')||document.querySelector('.popover[open],.fragment-settings[open],#coordinates-editor[open]')||!$('jobs-drawer').hidden||document.body.classList.contains('focus-view'))return;
      this.clear();event.preventDefault();
    },true);
  }
  kind(){return $('scan-coordinate').value;}
  arity(){return COORDINATES[this.kind()].arity;}
  indices(){return ['a','b','c','d'].slice(0,this.arity()).map(id=>$(`scan-${id}`).value===''?NaN:Number($(`scan-${id}`).value)-1);}
  complete(){const ids=this.indices(),m=this.getMolecule();return !!m&&ids.every(i=>Number.isInteger(i)&&m.atoms[i])&&new Set(ids).size===this.arity();}
  setIndices(ids){['a','b','c','d'].forEach((id,i)=>$(`scan-${id}`).value=Number.isInteger(ids[i])?ids[i]+1:'');}
  setKind(kind) {
    if(this.isBusy()||this.kind()===kind)return;
    $('scan-coordinate').value=kind;$('scan-unit').value=kind==='bond'?'angstrom':'degree';
    this.setIndices([]);$('scan-moving').value='';this.fragmentError='';
    $('scan-start').value=kind==='bond'?.7:kind==='angle'?80:0;$('scan-end').value=kind==='bond'?1.6:kind==='angle'?140:360;
    this.changed();
  }
  clear(){if(this.isBusy())return;this.setIndices([]);$('scan-moving').value='';this.fragmentError='';this.changed();}
  pick(index) {
    if(this.isBusy()||$('mode').value!=='scan')return;
    let ids=this.indices().filter(Number.isInteger);
    if(ids.length>=this.arity())ids=[];
    if(ids.includes(index))return;
    ids.push(index);this.setIndices(ids);
    if(this.complete())this.autoMoving();
    this.changed();
  }
  pickBond(a,b){if(this.kind()!=='bond'||this.isBusy()||$('mode').value!=='scan')return;this.setIndices([a,b]);this.autoMoving();this.resetRange();}
  autoMoving() {
    if(!this.complete())return;
    try {
      const m=this.getMolecule(),ids=this.indices(),kind=this.kind();
      const fixed=kind==='bond'?ids[0]:ids[1],movingSide=kind==='bond'?ids[1]:ids[2];
      const fragment=connectedFragment(m,fixed,movingSide);
      const forbidden=kind==='dihedral'?ids.slice(0,2):ids.slice(0,-1);
      if(forbidden.some(i=>fragment.includes(i))||!fragment.includes(ids.at(-1)))throw new Error('Choose moving atoms explicitly for this coordinate.');
      $('scan-moving').value=fragment.map(i=>i+1).join(', ');this.fragmentError='';
    }catch(error){$('scan-moving').value='';this.fragmentError=error.message;$('scan-picking-hint').textContent=error.message;}
  }
  current(){return coordinateValue(this.getMolecule(),this.kind(),this.indices(),$('scan-unit').value);}
  resetRange() {if(this.isBusy()||!this.complete())return;this.preset(0);}
  preset(index) {
    if(this.isBusy())return;
    const type=this.kind(),preset=COORDINATES[type].presets[index];let [,start,end,bracket]=preset;
    if(bracket!=null) {
      if(!this.complete())return;
      try{let delta=bracket;if(type==='bond'&&$('scan-unit').value==='bohr')delta/=ANGSTROM_PER_BOHR;const current=this.current();start=current-delta;end=current+delta;}catch{return;}
      if(type==='angle'){start=Math.max(.1,start);end=Math.min(180,end);}else start=Math.max(.01,start);
    }else if(type==='bond'&&$('scan-unit').value==='bohr'){start/=ANGSTROM_PER_BOHR;end/=ANGSTROM_PER_BOHR;}
    $('scan-start').value=start.toFixed(4);$('scan-end').value=end.toFixed(4);this.changed();
  }
  sync() {
    const type=this.kind(),meta=COORDINATES[type],m=this.getMolecule(),ids=this.indices();
    document.querySelectorAll('[data-coordinate]').forEach(button=>{button.setAttribute('aria-pressed',String(button.dataset.coordinate===type));button.disabled=this.isBusy();});
    for(const id of ['scan-reverse','scan-clear','scan-reset','scan-preset'])$(id).disabled=this.isBusy();
    document.querySelectorAll('.atom-token').forEach((token,i)=>{token.hidden=i>=meta.arity;token.querySelector('.atom-symbol').textContent=m?.atoms[ids[i]]?.symbol??'·';});
    document.querySelectorAll('#scan-unit option').forEach(option=>option.hidden=type==='bond'?option.value==='degree':option.value!=='degree');
    if(type==='bond'&&$('scan-unit').value==='degree')$('scan-unit').value='angstrom';
    if(type!=='bond')$('scan-unit').value='degree';
    if(this.lastKind!==type){$('scan-preset').replaceChildren(new Option('Scientific presets…',''),...meta.presets.map((p,i)=>new Option(p[0],i)));this.lastKind=type;}
    const chosen=ids.filter(i=>m?.atoms[i]);
    const label=i=>`${m.atoms[i].symbol}${i+1}`;
    let current='';
    if(this.complete())try{current=`${meta.symbol}₀ = ${this.current().toFixed(type==='bond'?4:2).replace('-', '−')} ${type==='bond'?($('scan-unit').value==='bohr'?'bohr':'Å'):'°'}`;}catch(error){current=error.message;}
    $('scan-current').textContent=current||`${meta.arity} atoms required`;
    $('selection-info').hidden=$('mode').value!=='scan'||!this.complete();
    $('selection-info').textContent=`${chosen.map(label).join('–')} · ${current}`;
    const next=chosen.length+1;
    $('scan-picking-hint').textContent=($('scan-moving').value.trim()?'':this.fragmentError)|| (this.complete()?'Coordinate ready. Click another atom to start a new selection.':`Click to assign Atom ${Math.min(next,meta.arity)} of ${meta.arity}.`);
    const start=Number($('scan-start').value),end=Number($('scan-end').value),count=Number($('scan-count').value),delta=(end-start)/(count-1);
    $('scan-delta').textContent=Number.isFinite(delta)&&count>=2?`Δ = ${delta>=0?'+':'−'}${Math.abs(delta).toFixed(4)} ${type==='bond'?($('scan-unit').value==='bohr'?'bohr':'Å'):'°'}/step`:'Δ = …';
    return chosen;
  }
}
