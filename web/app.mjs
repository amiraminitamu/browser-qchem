// SPDX-License-Identifier: AGPL-3.0-only
import {parseGeometry,packBasis,basisSummary,toXYZ,WATER_XYZ,ANGSTROM_PER_BOHR} from './molecule.mjs';
import {basisLabel} from './basis.mjs';
import {serializable} from './engine.mjs';
import {MoleculeViewer} from './viewer.mjs';
import {distance,connectedFragment,scanSpec,scanCSV,scanXYZ,HARTREE_TO_KCAL} from './scan.mjs';
import {plot} from './plot.mjs';

const $=id=>document.getElementById(id);
const EXAMPLES={water:WATER_XYZ,h2:'2\nHydrogen / angstrom\nH 0 0 0\nH 0 0 0.74',
  methane:'5\nMethane / angstrom\nC 0 0 0\nH .63 .63 .63\nH -.63 -.63 .63\nH -.63 .63 -.63\nH .63 -.63 -.63',
  ammonia:'4\nAmmonia / angstrom\nN 0 0 .1\nH .93 0 -.25\nH -.465 .8054 -.25\nH -.465 -.8054 -.25'};
let ready=false,busy=false,valid=false,jobId=0,jobMode='rhf',molecule=null,current=null;
let scan=null,selectedPoint=-1,history=[],previewTimer=null;
let pair=[0,1],selectionExplicit=false;
const viewer=new MoleculeViewer($('molecule-viewer'),selectAtom);
const worker=new Worker(new URL('./worker.mjs',import.meta.url),{type:'module'});
$('geometry').value=WATER_XYZ;

function setBusy(value) {
  busy=value;
  $('controls').querySelectorAll('input,select,textarea,button').forEach(element=>element.disabled=value);
  $('cancel').disabled=!value; $('cancel').textContent='Stop calculation';
  $('calculate').disabled=value||!ready||!valid;
  $('calculate').hidden=value;
  $('download').disabled=value||!current;
  $('xyz-file').disabled=value;
  $('geometry-download').disabled=!molecule;
  const mode=$('mode').value;
  document.querySelector('#calculation-section .method-chip').textContent=mode==='integrals'?'AO':'RHF';
  document.querySelector('.spin-info strong').textContent=mode==='integrals'?'Any electron count':'Closed-shell singlet';
  $('run-label').textContent=mode==='scan'?'Run bond scan':'Run calculation';
  $('mobile-run').disabled=!value&&(!ready||!valid);
  $('mobile-run').textContent=value?'Stop calculation':mode==='scan'?'Run bond scan':'Run calculation';
  $('scan-controls').querySelectorAll('input,select,button').forEach(element=>element.disabled=value||mode!=='scan');
  for(const id of ['charge','maxiter','diis'])$(id).disabled=value||mode==='integrals';
}
function badge(text,kind='neutral') { $('result-badge').textContent=text;$('result-badge').className=`status-badge ${kind}`; }
function setTab(name) {
  document.querySelectorAll('[data-tab]').forEach(button=>{
    button.setAttribute('aria-selected',String(button.dataset.tab===name));
    button.tabIndex=button.dataset.tab===name?0:-1;
  });
  document.querySelectorAll('[data-panel]').forEach(panel=>panel.hidden=panel.dataset.panel!==name);
  if(name==='scf')drawSCF();
  if(name==='pes')renderScan();
}
document.querySelectorAll('[data-tab]').forEach(button=>button.addEventListener('click',()=>setTab(button.dataset.tab)));
document.querySelectorAll('[data-tab]').forEach(button=>button.addEventListener('keydown',event=>{
  const tabs=[...document.querySelectorAll('[data-tab]')],index=tabs.indexOf(button);
  const target=event.key==='ArrowRight'?tabs[(index+1)%tabs.length]:event.key==='ArrowLeft'?tabs[(index+tabs.length-1)%tabs.length]:event.key==='Home'?tabs[0]:event.key==='End'?tabs.at(-1):null;
  if(target){event.preventDefault();setTab(target.dataset.tab);target.focus();}
}));
function formula(atoms) {
  const counts={};for(const atom of atoms) counts[atom.symbol]=(counts[atom.symbol]??0)+1;
  const symbols=Object.keys(counts).sort();
  if(counts.C) symbols.sort((a,b)=>(a==='C'?-2:a==='H'?-1:0)-(b==='C'?-2:b==='H'?-1:0)||a.localeCompare(b));
  return symbols.map(symbol=>symbol+(counts[symbol]>1?String(counts[symbol]).replace(/\d/g,d=>'₀₁₂₃₄₅₆₇₈₉'[d]):'')).join('');
}
function readMoving() {
  const text=$('scan-moving').value.trim();
  if(!text) return [];
  return text.split(/[,\s]+/).map(value=>/^\d+$/.test(value)?Number(value)-1:NaN);
}
function selection() {
  if(!molecule) return;
  const a=Number($('scan-a').value)-1,b=Number($('scan-b').value)-1;
  pair=[a,b];
  if(a===b||!molecule.atoms[a]||!molecule.atoms[b]) {
    $('selection-info').textContent='Select two different atoms to measure a bond.';viewer.select([]);return;
  }
  const label=index=>`${molecule.atoms[index].symbol}${index+1}`;
  $('selection-info').textContent=`${label(a)} - ${label(b)} · ${(distance(molecule,a,b)*ANGSTROM_PER_BOHR).toFixed(4)} Å`;
  const highlight=selectionExplicit||$('mode').value==='scan';
  viewer.select(highlight?pair:[],highlight?readMoving().filter(i=>molecule.atoms[i]):[]);
}
function selectAtom(index) {
  if(busy||!molecule||molecule.atoms.length<2) return;
  selectionExplicit=true;
  pair=pair.length===1?[pair[0],index]:[index];
  if(pair.length===1){viewer.select(pair);$('selection-info').textContent=`Atom ${index+1} selected. Choose a second atom.`;return;}
  if(pair[0]===pair[1]){pair=[index];return;}
  $('scan-a').value=pair[0]+1;$('scan-b').value=pair[1]+1;$('scan-moving').value=pair[1]+1;
  const d=distance(molecule,...pair)*ANGSTROM_PER_BOHR;
  $('scan-unit').value='angstrom';$('scan-start').value=(d*.75).toFixed(3);$('scan-end').value=(d*1.65).toFixed(3);
  selection();
}
function preview({reset=false}={}) {
  clearTimeout(previewTimer);
  try {
    molecule=parseGeometry($('geometry').value,$('unit').value,$('basis').value);
    const info=basisSummary(molecule);
    valid=info.fits&&($('mode').value!=='scan'||molecule.atoms.length>1);
    $('basis-info').classList.toggle('over-limit',!info.fits);
    $('ao-count').textContent=`${info.nao} atomic orbitals${!info.fits?' / limit 40':''}`;
    $('eri-memory').textContent=`${(info.eriBytes/1e6).toFixed(2)} MB ERIs`;
    $('input-message').textContent=!info.fits?'This basis exceeds the 40 AO limit. Select a smaller basis or molecule.':'';
    $('formula').textContent=formula(molecule.atoms);
    const title={water:'Water',h2:'Hydrogen',methane:'Methane',ammonia:'Ammonia'}[$('preset').value]??'Custom molecule';
    $('header-molecule').textContent=title;
    $('workspace-title').replaceChildren(document.createTextNode(title));
    const dot=document.createElement('span');dot.className='title-dot';dot.textContent='.';$('workspace-title').append(dot);
    $('workspace-description').textContent=$('mode').value==='scan'?'Map the energy along a molecular coordinate.':'Explore structure, energy, and electronic properties.';
    $('viewer-meta').textContent=`${molecule.atoms.length} atoms · ${basisLabel(molecule.basis)}`;
    if(!current){$('metric-basis').textContent=basisLabel(molecule.basis);$('metric-aos').textContent=info.nao;$('metric-cycles').textContent='...';}
    viewer.setMolecule(molecule,{reset});
    if(pair.some(i=>!molecule.atoms[i])) {$('scan-a').value=1;$('scan-b').value=Math.min(2,molecule.atoms.length);$('scan-moving').value=Math.min(2,molecule.atoms.length);}
    selection();renderData();
  } catch(error) {
    valid=false;$('input-message').textContent=error.message;
    $('coordinates-editor').open=true;
    $('viewer-meta').textContent='Showing last valid geometry';
  }
  $('scan-controls').hidden=$('mode').value!=='scan';
  if(!busy) setBusy(false);
}
function invalidate() {
  if(busy) return;
  current=null;$('download').disabled=true;badge('READY');$('energy-value').textContent='...';
  $('job-progress').hidden=true;
  $('status').textContent='Settings changed. Ready to calculate.';
  $('summary').textContent='Run a calculation for the current settings.';
  $('scf-history').replaceChildren();$('orbital-table').replaceChildren();$('matrix').replaceChildren();
  $('frontier-info').replaceChildren();$('ao-order').textContent='';$('eri-selectors').replaceChildren();$('eri-value').textContent='Run a single point to inspect integrals.';
  history=[];plot($('scf-chart'),[]);
  previewTimer=setTimeout(()=>preview(),120);
}
function table(container,headers,rows) {
  const t=document.createElement('table'),h=t.createTHead().insertRow(),body=t.createTBody();
  for(const text of headers){const th=document.createElement('th');th.textContent=text;h.append(th);}
  rows.forEach(values=>{const row=body.insertRow();values.forEach(value=>row.insertCell().textContent=value);});
  container.replaceChildren(t);return body;
}
function startHistory() { table($('scf-history'),['Cycle','Total E / Eh','delta E / Eh','Density RMS','Residual RMS'],[]); }
function appendIteration(entry) {
  const body=$('scf-history').querySelector('tbody')??(startHistory(),$('scf-history').querySelector('tbody'));
  const row=body.insertRow();row.dataset.iteration=entry.iteration;
  [entry.iteration,entry.energy.toFixed(12),entry.deltaEnergy.toExponential(3),entry.densityRMS.toExponential(3),entry.residualRMS.toExponential(3)].forEach(value=>row.insertCell().textContent=value);
  $('scf-history').scrollTop=$('scf-history').scrollHeight;
}
function drawSCF() {
  plot($('scf-chart'),history.map(entry=>({x:entry.iteration,y:entry.energy,label:`Cycle ${entry.iteration}: ${entry.energy.toFixed(12)} Eh`})),{
    xLabel:'SCF iteration',yLabel:'Energy / Eh',onSelect:index=>{
      $('iteration-details').open=true;
      $('scf-history').querySelectorAll('tbody tr').forEach((row,i)=>row.classList.toggle('selected',i===index));
      const row=$('scf-history').querySelectorAll('tbody tr')[index];
      if(row) $('scf-history').scrollTop=row.offsetTop-$('scf-history').offsetTop;
    }});
}
function renderOrbitals(rhf) {
  if(!rhf){$('orbital-table').textContent='Run RHF to view orbitals.';$('frontier-info').replaceChildren();return;}
  table($('orbital-table'),['Orbital','Energy / Eh','Occupation'],Array.from(rhf.orbitalEnergies,(e,i)=>[i+1,e.toFixed(10),rhf.occupations[i]]));
  const nocc=rhf.electrons/2,eps=rhf.orbitalEnergies;
  $('frontier-info').replaceChildren();
  for(const [label,value] of [['HOMO',eps[nocc-1]],['LUMO',eps[nocc]],['HOMO-LUMO gap',eps[nocc]===undefined?undefined:eps[nocc]-eps[nocc-1]]]) {
    const box=document.createElement('div');box.textContent=label;const strong=document.createElement('strong');strong.textContent=value===undefined?'Not in basis':`${value.toFixed(6)} Eh`;box.append(strong);$('frontier-info').append(box);
  }
}
function renderMatrix() {
  if(!current){$('matrix').textContent='Run a calculation to inspect its matrices.';return;}
  const {nao:n,aoLabels}=current,kind=$('matrix-kind').value;
  const data=kind==='core'?current.kinetic.map((value,i)=>value+current.nuclear[i]):['density','fock'].includes(kind)?current.rhf?.[kind]:current[kind];
  if(!data){$('matrix').textContent='Run RHF to inspect this matrix.';return;}
  table($('matrix'),['AO',...Array.from({length:n},(_,i)=>i+1)],Array.from({length:n},(_,i)=>[`${i+1}: ${aoLabels[i]}`,...Array.from({length:n},(_,j)=>data[i*n+j].toFixed(8))]));
}
function renderERI() {
  if(!current?.eri){$('eri-value').textContent='Run a single-point calculation at this geometry to inspect ERIs.';return;}
  const inputs=[...$('eri-selectors').querySelectorAll('input')];
  if(inputs.some(input=>!input.checkValidity()||!Number.isInteger(Number(input.value)))){$('eri-value').textContent='Choose four valid integer AO indices.';return;}
  const [i,j,k,l]=inputs.map(input=>Number(input.value)-1),n=current.nao;
  $('eri-value').textContent=`(${i+1} ${j+1}|${k+1} ${l+1}) = ${current.eri[((i*n+j)*n+k)*n+l].toFixed(14)} Eh`;
}
function renderData() {
  if(!molecule) return;
  $('geometry-output').textContent=toXYZ(molecule,'angstrom',`${basisLabel(molecule.basis)} / current geometry`);
  const info=basisSummary(molecule);$('data-meta').replaceChildren();
  for(const text of [`${molecule.atoms.length} atoms`,basisLabel(molecule.basis),`${info.nao} AOs`,'Coordinates exported in angstrom']){const chip=document.createElement('span');chip.textContent=text;$('data-meta').append(chip);}
}
function showResult(result,{pointIndex=null}={}) {
  current=result;const r=result.rhf,n=result.nao;
  $('metric-basis').textContent=basisLabel(result.molecule.basis);$('metric-aos').textContent=n;$('metric-cycles').textContent=r?.iterations??'N/A';
  if(r&&!r.converged)document.querySelector('.summary-details').open=true;
  $('result-context').textContent=`${r?'RHF':'INTEGRALS'} / ${basisLabel(result.molecule.basis)}`;
  $('results-title').textContent=pointIndex===null?'Calculation details':`Scan point ${pointIndex+1} / ${scan.spec.count}`;
  if(r){
    badge(r.converged?'CONVERGED':'NOT CONVERGED',r.converged?'success':'warning');
    $('energy-label').textContent=r.converged?'TOTAL RHF ENERGY':'UNCONVERGED ENERGY';$('energy-value').textContent=r.energy.toFixed(10);$('energy-unit').textContent='Hartree';
    $('status').textContent=r.converged?`CONVERGED: RHF/${basisLabel(result.molecule.basis)}.`:'NOT CONVERGED: maximum iterations reached.';
    $('summary').textContent=[r.converged?`Total RHF energy  ${r.energy.toFixed(12)} Eh`:`Unconverged energy ${r.energy.toFixed(12)} Eh`,
      `Electronic energy ${r.electronicEnergy.toFixed(10)}`,`Nuclear repulsion ${result.nuclearRepulsion.toFixed(10)}`,
      `Electrons / charge ${r.electrons} / ${r.charge}`,`Cycles / DIIS      ${r.iterations} / ${r.diis?'on':'off'}`,
      `Basis / AOs        ${basisLabel(result.molecule.basis)} / ${n}`,
      `Density RMS        ${r.history.at(-1).densityRMS.toExponential(2)}`,`Residual RMS       ${r.history.at(-1).residualRMS.toExponential(2)}`].join('\n');
  }else{
    badge(result.checks.passed?'INTEGRALS PASS':'CHECK FAILED',result.checks.passed?'success':'error');$('energy-label').textContent='ATOMIC ORBITALS';$('energy-value').textContent=n;$('energy-unit').textContent=basisLabel(result.molecule.basis);
    $('status').textContent=result.checks.passed?'PASS: normalization and symmetry checks.':'FAIL: integral checks.';
    $('summary').textContent=[`Atoms / shells   ${result.molecule.atoms.length} / ${result.nshell}`,`One-electron arrays 3 x ${n} x ${n}`,`ERI tensor ${n} x ${n} x ${n} x ${n} (${n**4} values)`,
      `Nuclear repulsion ${result.nuclearRepulsion.toFixed(12)} Eh`,`Max |Sii-1| ${result.checks.diagonalError.toExponential(2)}`,`ERI symmetry ${result.checks.eriSymmetryError.toExponential(2)}`].join('\n');
  }
  history=r?.history??[];startHistory();history.forEach(appendIteration);drawSCF();renderOrbitals(r);renderMatrix();
  $('ao-order').textContent=result.aoLabels.map((label,i)=>`${i+1}: ${label}`).join('\n');$('eri-selectors').replaceChildren();
  if(result.eri) for(const name of ['i','j','k','l']) {
    const label=document.createElement('label');label.textContent=name;const input=document.createElement('input');
    Object.assign(input,{type:'number',min:'1',max:String(n),step:'1',value:'1',required:true});input.setAttribute('aria-label',`ERI index ${name}`);
    input.addEventListener('input',renderERI);label.append(input);$('eri-selectors').append(label);
  }
  renderERI();$('download').disabled=busy||!current;
}
function scanConfig() {
  return {a:Number($('scan-a').value)-1,b:Number($('scan-b').value)-1,moving:readMoving(),start:Number($('scan-start').value),end:Number($('scan-end').value),count:Number($('scan-count').value),unit:$('scan-unit').value};
}
function options() {return {charge:Number($('charge').value),maxIterations:Number($('maxiter').value),diis:$('diis').checked};}
function run(event) {
  event?.preventDefault();if(busy||!ready) return;
  preview();if(!valid) return;
  if($('mode').value!=='integrals'){
    if(!$('maxiter').checkValidity())document.querySelector('.advanced-settings').open=true;
    if(!$('charge').reportValidity()||!$('maxiter').reportValidity())return;
  }
  const mode=$('mode').value;
  let spec;
  try {if(mode==='scan')spec=scanSpec(molecule,scanConfig());}
  catch(error){$('input-message').textContent=error.message;return;}
  jobId++;jobMode=mode;current=null;history=[];startHistory();drawSCF();renderOrbitals(null);renderMatrix();$('download').disabled=true;
  $('energy-value').textContent='...';badge('RUNNING');$('status').textContent='Computing molecular integrals...';$('summary').textContent='Running locally in your browser.';
  document.querySelector('.summary-details').open=false;
  $('result-context').textContent=`${mode==='integrals'?'INTEGRALS':'RHF'} / ${basisLabel(molecule.basis)}`;
  $('job-progress').hidden=false;$('progress-fill').style.width='0%';$('progress-text').textContent=mode==='scan'?`0 / ${spec.count} scan points`:'Starting SCF';
  if(mode==='scan') {scan={spec,source:structuredClone(molecule),options:options(),points:[],cancelled:false};selectedPoint=-1;renderScan();setTab('pes');}
  else setTab(mode==='integrals'?'matrices':'scf');
  setBusy(true);
  if(matchMedia('(max-width:850px)').matches)toggleSetup(false);
  worker.postMessage({id:jobId,geometry:$('geometry').value,unit:$('unit').value,basis:$('basis').value,mode,options:options(),scan:spec});
}
function selectPoint(index) {
  const point=scan?.points[index];if(!point)return;selectedPoint=index;
  $('geometry').value=toXYZ(point.molecule,'angstrom',`Scan point ${index+1}`);$('unit').value='angstrom';$('basis').value=scan.spec.basis;
  $('preset').value='custom';
  $('charge').value=scan.options.charge;$('maxiter').value=scan.options.maxIterations;$('diis').checked=scan.options.diis;
  preview();
  if(point.result)showResult(point.result,{pointIndex:index});
  else {
    current=null;history=[];startHistory();drawSCF();renderOrbitals(null);renderMatrix();
    $('ao-order').textContent='';$('eri-selectors').replaceChildren();renderERI();
    $('results-title').textContent=`Scan point ${index+1} / ${scan.spec.count}`;
    $('metric-cycles').textContent='...';document.querySelector('.summary-details').open=true;
    badge('POINT FAILED','error');$('energy-value').textContent='...';$('status').textContent=`Scan point ${index+1} failed.`;$('summary').textContent=point.error;$('download').disabled=true;
  }
  renderScan();
}
function renderScan() {
  const points=scan?.points??[];$('scan-count-badge').textContent=points.length;
  for(const id of ['scan-csv','scan-xyz','scan-json'])$(id).disabled=!points.length;
  if(!scan){plot($('pes-chart'),[],{empty:'Your potential energy curve will appear here.'});return;}
  const finite=points.filter(p=>Number.isFinite(p.energy)),converged=finite.filter(p=>p.converged);
  const baseline=Math.min(...(converged.length?converged:finite).map(p=>p.energy));
  const unit=$('pes-energy-unit').value,relative=unit!=='absolute-eh',factor=unit==='relative-kcal'?HARTREE_TO_KCAL:1;
  const data=points.map(p=>({x:p.coordinate,y:p.energy===null?null:(p.energy-(relative?baseline:0))*factor,failed:!p.converged,
    label:`Point ${p.index+1}: ${p.coordinate.toFixed(4)} ${scan.spec.unit}; ${p.energy??'no energy'} Eh; ${p.status}`}));
  plot($('pes-chart'),data,{xLabel:`Bond length / ${scan.spec.unit==='angstrom'?'Å':'bohr'}`,yLabel:unit==='relative-kcal'?'ΔE / kcal mol⁻¹':relative?'ΔE / Eh':'Energy / Eh',onSelect:selectPoint,selected:selectedPoint,empty:'Scan points will appear as calculations complete.'});
  $('pes-note').textContent=`${basisLabel(scan.spec.basis)} · unrelaxed RHF · ${points.length}/${scan.spec.count} points${scan.cancelled?' · stopped':''}. ${relative?(converged.length?'Relative to the lowest converged point.':'No converged energies yet.'):''} Failed points are not connected.`;
  const body=table($('scan-point-table'),['Point',`Distance / ${scan.spec.unit==='angstrom'?'Å':'bohr'}`,'Energy / Eh','Status','Cycles'],points.map(p=>[p.index+1,p.coordinate.toFixed(6),p.energy===null?'Not available':p.energy.toFixed(12),p.status,p.iterations]));
  [...body.rows].forEach((row,index)=>{row.className=`clickable ${selectedPoint===index?'selected':''}`;row.tabIndex=0;row.setAttribute('aria-label',`Inspect scan point ${index+1}`);row.addEventListener('click',()=>selectPoint(index));row.addEventListener('keydown',e=>{if(e.key==='Enter')selectPoint(index);});});
}
worker.onmessage=({data})=>{
  if(data.type==='ready') {ready=true;$('engine-state').textContent='Local engine ready';$('engine-dot').classList.add('ready');preview();setBusy(false);run();return;}
  if(data.type==='fatal') {ready=false;setBusy(false);badge('ENGINE ERROR','error');$('status').textContent='Engine unavailable.';$('summary').textContent=data.message;return;}
  if(data.id!==jobId)return;
  if(data.type==='iteration') {
    if(jobMode==='scan') {
      $('progress-text').textContent=`Point ${data.index+1} / ${scan.spec.count} · SCF cycle ${data.entry.iteration}`;
      if(!current){$('status').textContent=`SCF iteration ${data.entry.iteration} at scan point ${data.index+1}`;$('energy-value').textContent=data.entry.energy.toFixed(10);}
    }else{
      history.push(data.entry);appendIteration(data.entry);drawSCF();
      $('status').textContent=`SCF iteration ${data.entry.iteration}: E = ${data.entry.energy.toFixed(10)} Eh`;$('energy-value').textContent=data.entry.energy.toFixed(10);
      $('progress-text').textContent=`SCF cycle ${data.entry.iteration}`;
    }return;
  }
  if(data.type==='point') {
    const follow=selectedPoint<0||selectedPoint===scan.points.length-1;
    scan.points.push(data.point);scan.spec=data.spec;
    $('progress-fill').style.width=`${100*scan.points.length/scan.spec.count}%`;
    if(follow)selectPoint(scan.points.length-1);else renderScan();return;
  }
  if(data.type==='scan-complete') {
    scan.cancelled=data.summary.cancelled;scan.summary=data.summary;setBusy(false);renderScan();
    $('progress-text').textContent=`${scan.cancelled?'Stopped':'Complete'} · ${scan.points.length} / ${scan.spec.count} points`;
    if(!scan.points.length){badge('CANCELLED','warning');$('status').textContent='Scan cancelled before a point completed.';}
    return;
  }
  setBusy(false);$('job-progress').hidden=true;
  if(data.type==='result'){showResult(data.result);return;}
  current=null;$('download').disabled=true;$('energy-value').textContent='...';$('metric-cycles').textContent='...';
  document.querySelector('.summary-details').open=data.type!=='cancelled';
  badge(data.type==='cancelled'?'CANCELLED':'INPUT ERROR',data.type==='cancelled'?'warning':'error');
  $('status').textContent=data.type==='cancelled'?'Calculation cancelled.':'Input or calculation error.';$('summary').textContent=data.message??'You can adjust the settings and run again.';
};
worker.onerror=event=>{ready=false;setBusy(false);badge('WORKER ERROR','error');$('status').textContent='The calculation worker failed.';$('summary').textContent=`${event.message}\nReload the page to restart.`;};
function save(name,content,type='text/plain') {
  const url=URL.createObjectURL(new Blob([content],{type})),a=document.createElement('a');a.href=url;a.download=name;a.click();setTimeout(()=>URL.revokeObjectURL(url),1000);
}
$('controls').addEventListener('submit',run);
$('controls').addEventListener('invalid',event=>{
  const details=event.target.closest('details');
  if(details)details.open=true;
},true);
$('cancel').addEventListener('click',()=>{worker.postMessage({type:'cancel',id:jobId});$('cancel').disabled=true;$('cancel').textContent='Stopping...';$('mobile-run').disabled=true;$('mobile-run').textContent='Stopping...';$('progress-text').textContent='Stopping after the current integral or SCF step...';});
$('controls').addEventListener('input',event=>{
  if(busy||event.target.id==='xyz-file')return;
  if(event.target.id==='geometry')$('preset').value='custom';
  invalidate();if(event.target.id.startsWith('scan-'))selection();
});
$('mode').addEventListener('change',()=>{
  if(busy)return;
  preview();
  if($('mode').value==='scan')$('scan-controls').scrollIntoView({block:'nearest',behavior:matchMedia('(prefers-reduced-motion:reduce)').matches?'instant':'smooth'});
});
$('preset').addEventListener('change',()=>{if(busy)return;selectionExplicit=false;$('geometry').value=EXAMPLES[$('preset').value];$('unit').value='angstrom';$('charge').value='0';invalidate();preview({reset:true});});
$('water').addEventListener('click',()=>{if(busy)return;$('preset').value='water';$('geometry').value=WATER_XYZ;$('unit').value='angstrom';$('charge').value='0';invalidate();preview({reset:true});run();});
$('xyz-file').addEventListener('change',async event=>{
  const file=event.target.files[0];if(!file)return;
  if(file.size>1000000){$('input-message').textContent='Use an XYZ file smaller than 1 MB.';return;}
  const text=await file.text();if(busy)return;
  $('geometry').value=text;$('preset').value='custom';$('unit').value='angstrom';invalidate();preview({reset:true});event.target.value='';
});
$('auto-fragment').addEventListener('click',()=>{try{$('scan-moving').value=connectedFragment(molecule,Number($('scan-a').value)-1,Number($('scan-b').value)-1).map(i=>i+1).join(', ');selection();$('input-message').textContent='Moving fragment highlighted. Check the selected atom indices before running.';}catch(e){$('input-message').textContent=e.message;}});
$('atom-labels').addEventListener('change',()=>viewer.setLabels($('atom-labels').checked));
$('fit-view').addEventListener('click',()=>viewer.reset());
$('viewer-png').addEventListener('click',()=>{const uri=viewer.png();if(!uri){$('input-message').textContent='Image export needs a working WebGL viewer.';return;}const a=document.createElement('a');a.href=uri;a.download='molecule.png';a.click();});
$('geometry-download').addEventListener('click',()=>{if(molecule)save('geometry.xyz',toXYZ(molecule,'angstrom',`${basisLabel(molecule.basis)} / exported geometry`));});
$('matrix-kind').addEventListener('change',renderMatrix);
$('pes-energy-unit').addEventListener('change',renderScan);
$('download').addEventListener('click',()=>{if(current)save('calculation.json',JSON.stringify(serializable(current),null,2),'application/json');});
$('scan-csv').addEventListener('click',()=>save('pes-scan.csv',scanCSV(scan),'text/csv'));
$('scan-xyz').addEventListener('click',()=>save('pes-scan.xyz',scanXYZ(scan)));
$('scan-json').addEventListener('click',()=>save('pes-scan.json',JSON.stringify({...scan,points:scan.points.map(p=>({...p,...(p.result?{result:serializable(p.result)}:{})}))},null,2),'application/json'));
// The setup drawer and viewer focus mode share the existing calculation state.
function toggleSetup(open) {
  const mobile=matchMedia('(max-width:850px)').matches;
  document.body.classList.toggle('setup-open',mobile&&open);
  $('setup-backdrop').hidden=!mobile||!open;
  $('mobile-setup-toggle').setAttribute('aria-expanded',String(mobile&&open));
  $('setup-panel').setAttribute('aria-hidden',String(mobile&&!open));
  if(mobile&&open)$('preset').focus();
  if(mobile&&!open&&$('setup-panel').contains(document.activeElement))$('mobile-setup-toggle').focus();
}
function focusViewer(enabled) {
  document.body.classList.toggle('focus-view',enabled);
  $('focus-view').setAttribute('aria-pressed',String(enabled));
  $('focus-view').title=enabled?'Exit expanded viewer':'Expand molecule viewer';
}
function activateTool(id) {
  document.querySelectorAll('.rail-button').forEach(button=>{
    button.classList.toggle('active',button.id===id);
    button.setAttribute('aria-pressed',String(button.id===id));
  });
}
$('mobile-setup-toggle').addEventListener('click',()=>toggleSetup(!document.body.classList.contains('setup-open')));
$('mobile-edit').addEventListener('click',()=>toggleSetup(true));
$('close-setup').addEventListener('click',()=>toggleSetup(false));
$('setup-backdrop').addEventListener('click',()=>toggleSetup(false));
$('mobile-run').addEventListener('click',()=>busy?$('cancel').click():run());
$('focus-view').addEventListener('click',()=>focusViewer(!document.body.classList.contains('focus-view')));
$('help-open').addEventListener('click',()=>$('workspace-guide').showModal());
for(const id of ['help-close','help-done'])$(id).addEventListener('click',()=>$('workspace-guide').close());
$('nav-molecule').addEventListener('click',()=>{
  activateTool('nav-molecule');$('coordinates-editor').open=true;toggleSetup(true);$('geometry').focus();
});
$('nav-setup').addEventListener('click',()=>{activateTool('nav-setup');toggleSetup(true);$('basis').focus();});
$('nav-scan').addEventListener('click',()=>{
  if(busy)return;
  activateTool('nav-scan');$('mode').value='scan';invalidate();preview();setTab('pes');toggleSetup(true);$('scan-a').focus();
});
$('nav-results').addEventListener('click',()=>{
  activateTool('nav-results');$('results-section').scrollIntoView({block:'start',behavior:matchMedia('(prefers-reduced-motion:reduce)').matches?'instant':'smooth'});
});
document.addEventListener('keydown',event=>{
  if(event.key==='Escape'){focusViewer(false);toggleSetup(false);}
  if((event.metaKey||event.ctrlKey)&&event.key==='Enter'&&!$('workspace-guide').open&&!busy){event.preventDefault();run();}
  if(event.key==='Tab'&&!$('workspace-guide').open){
    const root=document.body.classList.contains('focus-view')?document.querySelector('.viewer-card'):document.body.classList.contains('setup-open')?$('setup-panel'):null;
    if(root){
      const elements=[...root.querySelectorAll('button,input,select,textarea,[tabindex]')].filter(element=>!element.disabled&&element.tabIndex>=0&&element.offsetParent!==null);
      const first=elements[0],last=elements.at(-1);
      if(first&&(!root.contains(document.activeElement)||event.shiftKey&&document.activeElement===first||!event.shiftKey&&document.activeElement===last)){
        event.preventDefault();(event.shiftKey?last:first).focus();
      }
    }
  }
});
matchMedia('(max-width:850px)').addEventListener('change',()=>toggleSetup(false));
const chartWidths=new WeakMap();
new ResizeObserver(entries=>{
  for(const entry of entries){
    const width=Math.round(entry.contentRect.width);
    if(width>0&&chartWidths.get(entry.target)!==width){
      chartWidths.set(entry.target,width);
      if(entry.target.id==='scf-chart')drawSCF();else renderScan();
    }
  }
}).observe($('scf-chart'));
// Observe the second chart with the same width-only rule to avoid resize feedback.
new ResizeObserver(entries=>{
  for(const entry of entries){const width=Math.round(entry.contentRect.width);if(width>0&&chartWidths.get(entry.target)!==width){chartWidths.set(entry.target,width);renderScan();}}
}).observe($('pes-chart'));
toggleSetup(false);

preview({reset:true});renderScan();plot($('scf-chart'),[]);renderOrbitals(null);renderMatrix();setBusy(false);

