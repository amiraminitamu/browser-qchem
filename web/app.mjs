// SPDX-License-Identifier: AGPL-3.0-only
import {parseGeometry,packBasis,basisSummary,toXYZ,WATER_XYZ,ANGSTROM_PER_BOHR} from './molecule.mjs';
import {basisLabel} from './basis.mjs';
import {serializable} from './engine.mjs';
import {MoleculeViewer} from './viewer.mjs';
import {distance,scanSpec,scanCSV,scanXYZ} from './scan.mjs';
import {scfResult,totalEnergy,validateElectronicOptions} from './electronic.mjs';
import {PESSetup} from './pes-setup.mjs';
import {COORDINATES} from './coordinates.mjs';
import {exportPESFigure,pesFigureModel} from './pes-export.mjs';
import {mathText} from './format.mjs';
import {plot} from './plot.mjs';
import {BUILD_ID,WORKER_PROTOCOL} from './build-info.mjs';
import {initializeWorkbench,toggleSetup,TASKS} from './workbench.mjs';
import {captureRequest,assertResultMatches,presetMethod,restoreSettings,saveSettings,readSettings,SETTINGS_FIELDS} from './settings.mjs';

const $=id=>document.getElementById(id);
initializeWorkbench();
const EXAMPLES={water:WATER_XYZ,h2:'2\nHydrogen / angstrom\nH 0 0 0\nH 0 0 0.74',
  methane:'5\nMethane / angstrom\nC 0 0 0\nH .63 .63 .63\nH -.63 -.63 .63\nH -.63 .63 -.63\nH .63 -.63 -.63',
  ammonia:'4\nAmmonia / angstrom\nN 0 0 .1\nH .93 0 -.25\nH -.465 .8054 -.25\nH -.465 -.8054 -.25',
  h2plus:'2\nHydrogen molecular ion / angstrom\nH 0 0 0\nH 0 0 1.06',
  heh:'2\nNeutral helium hydride / angstrom\nHe 0 0 0\nH 0 0 1.0',
  hehplus:'2\nHelium hydride cation / angstrom\nHe 0 0 0\nH 0 0 .8'};
const EXAMPLE_SETTINGS={h2plus:{name:'Hydrogen ion',charge:1,multiplicity:2,method:'uhf'},heh:{name:'Helium hydride',charge:0,multiplicity:2,method:'uhf'},hehplus:{name:'Helium hydride ion',charge:1,multiplicity:1,method:'rhf'}};
let ready=false,busy=false,valid=false,jobId=0,jobMode='rhf',molecule=null,current=null;
let scan=null,selectedPoint=-1,history=[],previewTimer=null,activeRequest=null,draftRevision=0,scanRevision=0,methodExplicit=false;
let listedJobs=[],listedJobSeq=0,activeListedJob=0;
const viewer=new MoleculeViewer($('molecule-viewer'),selectAtom,(a,b)=>pesSetup.pickBond(a,b));
const pesSetup=new PESSetup({getMolecule:()=>molecule,isBusy:()=>busy,changed:()=>settingsChanged()});
const worker=new Worker(new URL('./worker.mjs',import.meta.url),{type:'module'});
$('controls').reset();
$('geometry').value=WATER_XYZ;
function storage() {try {return globalThis.localStorage;} catch {return null;}}
methodExplicit=restoreSettings($('controls'),storage());
function persistSettings() {saveSettings($('controls'),storage(),methodExplicit);}
globalThis.addEventListener('pagehide',persistSettings);

function setDisabled(element,value) {if(element.disabled!==value)element.disabled=value;}
function setBusy(value) {busy=value;syncControls();}
function syncControls() {
  const mode=$('mode').value,uhf=$('method').value==='uhf';
  // Never reset all controls before locking a subset. In particular, leave a
  // native select alone when its enabled state has not actually changed.
  $('controls').querySelectorAll('input,select,textarea,button').forEach(element=>{
    if(!SETTINGS_FIELDS.includes(element.id)&&!['xyz-file','water','auto-fragment'].includes(element.id))return;
    const inactive=element.closest('#scan-controls')&&(mode!=='scan'||(['scan-c','scan-d'].includes(element.id)&&Number(element.id==='scan-c'?3:4)>({bond:2,angle:3,dihedral:4})[$('scan-coordinate').value]));
    const scf=mode==='integrals'&&['charge','maxiter','diis','multiplicity','guess'].includes(element.id);
    const spin=!uhf&&['multiplicity','guess'].includes(element.id);
    setDisabled(element,!!(busy||inactive||scf||spin));
  });
  if(typeof pesSetup!=='undefined')pesSetup.sync();
  setDisabled($('cancel'),!busy);$('cancel').hidden=!busy;$('cancel').innerHTML='<svg class="icon" aria-hidden="true"><use href="#i-stop"/></svg> Stop';
  setDisabled($('calculate'),busy||!ready||!valid);$('calculate').hidden=busy;
  setDisabled($('download'),busy||!current);
  setDisabled($('geometry-download'),!molecule);
  document.querySelector('#solver-popover .method-chip').textContent=mode==='integrals'?'AO':$('method').value.toUpperCase();
  document.querySelector('.spin-info strong').textContent=mode==='integrals'?'Any electron count':uhf?`Multiplicity ${$('multiplicity').value}`:'Closed-shell singlet';
  $('run-label').textContent='Run';
  $('calculate').title=mode==='scan'?'Run PES scan (Cmd/Ctrl + Enter)':'Run calculation (Cmd/Ctrl + Enter)';
  $('uhf-controls').hidden=!uhf||mode==='integrals';
  $('method-hint').textContent=uhf?'Separate alpha and beta orbitals, with spin diagnostics.':$('method').value==='mp2'?'All-electron MP2 correlation on an RHF reference.':'Restricted Hartree-Fock for closed-shell singlets.';
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
  if(symbols.length===2&&counts.He&&counts.H) symbols.sort((a,b)=>a==='He'?-1:b==='He'?1:0);
  if(counts.C) symbols.sort((a,b)=>(a==='C'?-2:a==='H'?-1:0)-(b==='C'?-2:b==='H'?-1:0)||a.localeCompare(b));
  return symbols.map(symbol=>symbol+(counts[symbol]>1?String(counts[symbol]).replace(/\d/g,d=>'₀₁₂₃₄₅₆₇₈₉'[d]):'')).join('');
}
function chargeSuffix(charge) {
  const digits=String(Math.abs(charge)).replace(/\d/g,d=>'⁰¹²³⁴⁵⁶⁷⁸⁹'[d]);
  return (Math.abs(charge)===1?'':digits)+(charge>0?'⁺':'⁻');
}
function readMoving() {
  const text=$('scan-moving').value.trim();
  if(!text) return [];
  return text.split(/[,\s]+/).map(value=>/^\d+$/.test(value)?Number(value)-1:NaN);
}
function selection() {
  if(!molecule)return;
  const task=TASKS.find(t=>t.id===$('mode').value);
  viewer.setPickingMode(task.picking,$('scan-coordinate').value);
  const ids=pesSetup.sync();
  viewer.select(task.picking==='internal-coordinate'?ids:[],task.picking==='internal-coordinate'?readMoving().filter(i=>molecule.atoms[i]):[]);
}
function selectAtom(index) {pesSetup.pick(index);}
function preview({reset=false}={}) {
  clearTimeout(previewTimer);
  try {
    molecule=parseGeometry($('geometry').value,$('unit').value,$('basis').value);
    const info=basisSummary(molecule);
    valid=info.fits&&($('mode').value!=='scan'||molecule.atoms.length>1);
    $('basis-info').classList.toggle('over-limit',!info.fits);
    $('ao-count').textContent=`${info.nao} AOs${!info.fits?' / limit 40':''}`;
    $('eri-memory').textContent=`${(info.eriBytes/1e6).toFixed(2)} MB ERIs`;
    $('input-message').textContent=!info.fits?'This basis exceeds the 40 AO limit. Select a smaller basis or molecule.':'';
    $('formula').textContent=formula(molecule.atoms);
    const title=EXAMPLE_SETTINGS[$('preset').value]?.name??{water:'Water',h2:'Hydrogen',methane:'Methane',ammonia:'Ammonia'}[$('preset').value]??'Custom molecule';
    $('header-molecule').textContent=title;
    $('workspace-title').replaceChildren(document.createTextNode(title));
    const charge=Number($('charge').value),electrons=molecule.atoms.reduce((sum,a)=>sum+a.charge,0)-charge;
    const mult=$('method').value==='uhf'?Number($('multiplicity').value):1;
    const spinLabel=({1:'singlet',2:'doublet',3:'triplet',4:'quartet',5:'quintet'})[mult]??`multiplicity ${mult}`;
    const chargeLabel=Number.isSafeInteger(charge)?charge===0?'Neutral':`Charge ${charge>0?'+':''}${charge}`:'Invalid charge';
    $('molecule-state').textContent=mathText($('mode').value==='integrals'?chargeLabel:`${chargeLabel} · ${spinLabel}`);
    $('electron-info').textContent=Number.isSafeInteger(electrons)&&electrons>=0?`${electrons} e⁻`:'Check the molecular charge';
    $('formula').textContent=formula(molecule.atoms)+(Number.isSafeInteger(charge)&&charge!==0?chargeSuffix(charge):'');
    $('workspace-description').textContent=$('mode').value==='scan'?'Map the energy along a molecular coordinate.':'Explore structure, energy, and electronic properties.';
    $('viewer-meta').textContent=`${molecule.atoms.length} atoms · ${basisLabel(molecule.basis)}`;
    if(!current){$('metric-basis').textContent=basisLabel(molecule.basis);$('metric-aos').textContent=info.nao;$('metric-cycles').textContent='...';$('metric-electrons').textContent=Number.isSafeInteger(electrons)&&electrons>=0?electrons:'...';}
    viewer.setMolecule(molecule,{reset});
    selection();renderData();
  } catch(error) {
    valid=false;$('input-message').textContent=error.message;
    $('coordinates-editor').open=true;
    $('viewer-meta').textContent='Showing last valid geometry';
  }
  $('scan-controls').hidden=TASKS.find(t=>t.id===$('mode').value)?.context!=='scan-controls';
  document.body.dataset.task=$('mode').value;
  syncControls();
}
function invalidate() {
  if(busy) return;
  clearTimeout(previewTimer);++draftRevision;++jobId;activeRequest=null;
  selectedPoint=-1;
  const method=$('method').value.toUpperCase();
  $('energy-label').textContent=$('mode').value==='integrals'?'ATOMIC ORBITALS':`TOTAL ${method} ENERGY`;
  $('energy-unit').textContent='Eₕ';
  $('result-context').textContent=`${$('mode').value==='integrals'?'INTEGRALS':method} / ${basisLabel($('basis').value)}`;
  referenceContext(method==='MP2');
  current=null;$('method-details').hidden=true;$('download').disabled=true;badge('READY');$('energy-value').textContent=mathText('...');
  $('job-progress').hidden=true;
  $('status').textContent='Settings changed. Ready to calculate.';
  $('summary').textContent='Run a calculation for the current settings.';
  $('scf-history').replaceChildren();$('orbital-table').replaceChildren();$('matrix').replaceChildren();
  $('frontier-info').replaceChildren();
  document.dispatchEvent(new CustomEvent('lorentzian:orbital-select',{detail:null}));
  $('ao-order').textContent='';$('eri-selectors').replaceChildren();$('eri-value').textContent='Run a single point to inspect integrals.';
  history=[];plot($('scf-chart'),[]);
  previewTimer=setTimeout(()=>preview(),120);
}
function table(container,headers,rows) {
  const t=document.createElement('table'),h=t.createTHead().insertRow(),body=t.createTBody();
  for(const text of headers){const th=document.createElement('th');th.textContent=mathText(text);h.append(th);}
  rows.forEach(values=>{const row=body.insertRow();values.forEach(value=>row.insertCell().textContent=mathText(value));});
  container.replaceChildren(t);return body;
}
function startHistory() { table($('scf-history'),['Cycle','Total E / Eₕ','ΔE / Eₕ','Density RMS','Residual RMS'],[]); }
function appendIteration(entry) {
  const body=$('scf-history').querySelector('tbody')??(startHistory(),$('scf-history').querySelector('tbody'));
  const row=body.insertRow();row.dataset.iteration=entry.iteration;
  [entry.iteration,entry.energy.toFixed(12),entry.deltaEnergy.toExponential(3),entry.densityRMS.toExponential(3),entry.residualRMS.toExponential(3)].forEach(value=>row.insertCell().textContent=mathText(value));
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
function renderOrbitals(r) {
  $('frontier-info').replaceChildren();
  document.dispatchEvent(new CustomEvent('lorentzian:orbital-select',{detail:null}));
  if(!r){$('orbital-table').textContent='Run an electronic calculation to view orbitals.';return;}
  const uhf=r.method==='UHF';
  $('orbital-note').textContent=uhf?'Separate alpha and beta orbital energies in Hartree, occupation 1 or 0.':current?.mp2?'RHF reference orbital energies and occupations. These are not correlated MP2 orbitals.':'Spatial orbital energies in Hartree. Closed-shell occupation is 2 or 0.';
  const channels=uhf?[['Alpha',r.orbitalEnergiesAlpha,r.occupationsAlpha,r.nalpha],['Beta',r.orbitalEnergiesBeta,r.occupationsBeta,r.nbeta]]:[['',r.orbitalEnergies,r.occupations,r.electrons/2]];
  const rows=[];
  for(const [label,eps,occ,nocc] of channels) {
    eps.forEach((e,i)=>rows.push({spin:label,index:i,energy:e,values:[...(uhf?[label]:[]),i+1,occ[i],e.toFixed(10),(e*27.211386245988).toFixed(6)]}));
    for(const [name,value] of [['HOMO',eps[nocc-1]],['LUMO',eps[nocc]],['Gap',nocc>0&&eps[nocc]!==undefined?eps[nocc]-eps[nocc-1]:undefined]]) {
      const box=document.createElement('div');box.textContent=`${label} ${name}`.trim();const strong=document.createElement('strong');strong.textContent=value===undefined?'Not in basis':mathText(`${value.toFixed(6)} Eₕ`);box.append(strong);$('frontier-info').append(box);
    }
  }
  const body=table($('orbital-table'),[...(uhf?['Spin']:[]),'Orbital','Occupation','Energy / Eh','Energy / eV'],rows.map(row=>row.values));
  [...body.rows].forEach((row,index)=>{
    const orbital=rows[index];row.className='clickable';row.tabIndex=0;row.setAttribute('aria-selected','false');
    row.setAttribute('aria-label',`${orbital.spin} orbital ${orbital.index+1}, ${orbital.energy.toFixed(6)} Hartree`);
    function select(){
      [...body.rows].forEach(other=>{other.classList.toggle('selected',other===row);other.setAttribute('aria-selected',String(other===row));});
      document.dispatchEvent(new CustomEvent('lorentzian:orbital-select',{detail:orbital}));
    }
    row.addEventListener('click',select);row.addEventListener('keydown',event=>{if(event.key==='Enter'||event.key===' '){event.preventDefault();select();}});
  });
}
function renderMatrix() {
  if(!current){$('matrix').textContent='Run a calculation to inspect its matrices.';return;}
  const {nao:n,aoLabels}=current,kind=$('matrix-kind').value;
  const data=kind==='core'?current.kinetic.map((value,i)=>value+current.nuclear[i]):['density','fock','densityAlpha','densityBeta','spinDensity','fockAlpha','fockBeta'].includes(kind)?scfResult(current)?.[kind]:current[kind];
  if(!data){$('matrix').textContent='This matrix is not available for this method. For UHF, choose alpha or beta Fock.';return;}
  table($('matrix'),['AO',...Array.from({length:n},(_,i)=>i+1)],Array.from({length:n},(_,i)=>[`${i+1}: ${aoLabels[i]}`,...Array.from({length:n},(_,j)=>data[i*n+j].toFixed(8))]));
}
function renderERI() {
  if(!current?.eri){$('eri-value').textContent='Run a single-point calculation at this geometry to inspect ERIs.';return;}
  const inputs=[...$('eri-selectors').querySelectorAll('input')];
  if(inputs.some(input=>!input.checkValidity()||!Number.isInteger(Number(input.value)))){$('eri-value').textContent='Choose four valid integer AO indices.';return;}
  const [i,j,k,l]=inputs.map(input=>Number(input.value)-1),n=current.nao;
  $('eri-value').textContent=mathText(`(${i+1} ${j+1}|${k+1} ${l+1}) = ${current.eri[((i*n+j)*n+k)*n+l].toFixed(14)} Eₕ`);
}
function renderData() {
  if(!molecule)return;
  $('geometry-output').textContent=toXYZ(molecule,'angstrom',`${basisLabel(molecule.basis)} / current geometry`);
  const info=basisSummary(molecule);$('data-meta').replaceChildren();
  for(const text of [`${molecule.atoms.length} atoms`,basisLabel(molecule.basis),`${info.nao} AOs`]){const chip=document.createElement('span');chip.textContent=text;$('data-meta').append(chip);}
  const atoms=molecule.atoms,label=(atom,i)=>`${atom.symbol}${i+1}`;
  table($('coordinate-table'),['Atom','x / Å','y / Å','z / Å'],atoms.map((atom,i)=>[label(atom,i),...atom.xyz.map(x=>(x*ANGSTROM_PER_BOHR).toFixed(8))]));
  table($('distance-matrix'),['Atom',...atoms.map(label)],atoms.map((atom,i)=>[label(atom,i),...atoms.map((_,j)=>(distance(molecule,i,j)*ANGSTROM_PER_BOHR).toFixed(5))]));
}
function referenceContext(mp2) {
  document.querySelector('#panel-scf h3').textContent=mp2?'RHF reference convergence':'SCF convergence';
  document.querySelector('#panel-scf .panel-intro p').textContent=mp2?'The SCF curve belongs to the RHF reference; correlation is computed after convergence.':'Follow the total SCF energy as the density converges.';
}
function renderMethodDetails(result) {
  referenceContext(!!result.mp2);
  const box=$('method-details');box.replaceChildren();const u=result.uhf,m=result.mp2;
  $('method-note').textContent='SCF convergence does not establish wavefunction stability or the lowest-energy state.';
  box.hidden=!u&&!m;if(box.hidden)return;
  const values=u?[['Alpha / beta electrons',`${u.nalpha} / ${u.nbeta}`],['Multiplicity',u.multiplicity],['⟨S²⟩ / target',`${u.spinSquare.toFixed(6)} / ${u.targetSpinSquare.toFixed(6)}`],['Spin contamination',u.spinContamination.toFixed(6)]]:
    m.status==='complete'?[['RHF reference',`${m.referenceEnergy.toFixed(10)} Eh`],['MP2 correlation',`${m.correlationEnergy.toFixed(10)} Eh`],['Opposite spin',`${m.oppositeSpinEnergy.toFixed(10)} Eh`],['Same spin',`${m.sameSpinEnergy.toFixed(10)} Eh`]]:[['MP2','Not computed: RHF did not converge']];
  for(const [label,value] of values){const item=document.createElement('div');item.textContent=label;const strong=document.createElement('strong');strong.textContent=mathText(value);item.append(strong);box.append(item);}
  $('method-note').textContent=u?`Guess: ${u.mixedGuessApplied?'opposite alpha/beta HOMO-LUMO mixing':'core Hamiltonian'}. UHF is not generally an eigenstate of total spin; convergence does not guarantee the lowest or stable SCF solution.`:'All-electron conventional MP2. Orbitals, density, Fock matrices and SCF history belong to the RHF reference. MP2 can be unreliable near bond dissociation or small orbital gaps.';
}
function showResult(result,{pointIndex=null,skipRemember=false}={}) {
  current=result;const r=scfResult(result),n=result.nao,method=result.method??r?.method??'INTEGRALS',energy=totalEnergy(result);
  renderMethodDetails(result);
  $('metric-basis').textContent=basisLabel(result.molecule.basis);$('metric-aos').textContent=n;$('metric-cycles').textContent=r?.iterations??'N/A';$('metric-electrons').textContent=r?.electrons??'N/A';
  if(r&&!r.converged)document.querySelector('.summary-details').open=true;
  $('result-context').textContent=`${method} / ${basisLabel(result.molecule.basis)}`;
  $('results-title').textContent=pointIndex===null?'Electronic structure':`Scan point ${pointIndex+1} / ${scan.spec.count}`;
  if(r){
    badge(r.converged?'CONVERGED':'NOT CONVERGED',r.converged?'success':'warning');
    $('energy-label').textContent=r.converged?`TOTAL ${method} ENERGY`:'UNCONVERGED ENERGY';$('energy-value').textContent=mathText((energy??r.energy).toFixed(10));$('energy-unit').textContent='Eₕ';
    $('status').textContent=r.converged?`CONVERGED: ${method}/${basisLabel(result.molecule.basis)}.`:result.mp2?'NOT CONVERGED: RHF reference; MP2 not computed.':'NOT CONVERGED: maximum iterations reached.';
    $('summary').textContent=mathText([r.converged?`Total ${method} energy  ${(energy??r.energy).toFixed(12)} Eh`:`Unconverged energy ${r.energy.toFixed(12)} Eh`,
      `Electronic energy ${r.electronicEnergy.toFixed(10)}`,`Nuclear repulsion ${result.nuclearRepulsion.toFixed(10)}`,
      `Electrons / charge ${r.electrons} / ${r.charge}`,`Cycles / DIIS      ${r.iterations} / ${r.diis?'on':'off'}`,
      `Basis / AOs        ${basisLabel(result.molecule.basis)} / ${n}`,
      `Density RMS        ${r.history.at(-1).densityRMS.toExponential(2)}`,`Residual RMS       ${r.history.at(-1).residualRMS.toExponential(2)}`].join('\n'));
  }else{
    badge(result.checks.passed?'INTEGRALS PASS':'CHECK FAILED',result.checks.passed?'success':'error');$('energy-label').textContent='ATOMIC ORBITALS';$('energy-value').textContent=mathText(n);$('energy-unit').textContent=basisLabel(result.molecule.basis);
    $('status').textContent=result.checks.passed?'PASS: normalization and symmetry checks.':'FAIL: integral checks.';
    $('summary').textContent=mathText([`Atoms / shells   ${result.molecule.atoms.length} / ${result.nshell}`,`One-electron arrays 3 x ${n} x ${n}`,`ERI tensor ${n} x ${n} x ${n} x ${n} (${n**4} values)`,
      `Nuclear repulsion ${result.nuclearRepulsion.toFixed(12)} Eh`,`Max |Sii-1| ${result.checks.diagonalError.toExponential(2)}`,`ERI symmetry ${result.checks.eriSymmetryError.toExponential(2)}`].join('\n'));
  }
  history=r?.history??[];startHistory();history.forEach(appendIteration);drawSCF();renderOrbitals(r);renderMatrix();
  $('ao-order').textContent=result.aoLabels.map((label,i)=>`${i+1}: ${label}`).join('\n');$('eri-selectors').replaceChildren();
  if(result.eri) for(const name of ['i','j','k','l']) {
    const label=document.createElement('label');label.textContent=name;const input=document.createElement('input');
    Object.assign(input,{type:'number',min:'1',max:String(n),step:'1',value:'1',required:true});input.setAttribute('aria-label',`ERI index ${name}`);
    input.addEventListener('input',renderERI);label.append(input);$('eri-selectors').append(label);
  }
  renderERI();$('download').disabled=busy||!current;
  if(!skipRemember&&pointIndex===null) rememberCompletedJob(result);
}
function run(event) {
  event?.preventDefault();if(busy||!ready)return;
  clearTimeout(previewTimer);
  let request,spec,target;
  try {
    request=captureRequest($('controls'),{id:jobId+1,protocol:WORKER_PROTOCOL,buildId:BUILD_ID});
    preview();if(!valid)return;
    if(request.mode!=='integrals') {
      if(!$('maxiter').checkValidity())document.querySelector('.advanced-settings').open=true;
      if(!$('charge').reportValidity()||!$('maxiter').reportValidity())return;
    }
    target=parseGeometry(request.geometry,request.unit,request.basis);
    if(request.mode!=='integrals')validateElectronicOptions(target,request.options);
    if(request.mode==='scan')spec={...scanSpec(target,request.scan),...request.options};
  } catch(error) {
    $('input-message').textContent=error.message;badge('INPUT ERROR','error');
    $('status').textContent='Input or calculation error.';$('summary').textContent=error.message;
    document.querySelector('.summary-details').open=true;return;
  }
  persistSettings();
  $('method-details').hidden=true;referenceContext(request.options.method==='mp2');
  jobId=request.id;jobMode=request.mode;activeRequest=request;current=null;history=[];
  startHistory();drawSCF();renderOrbitals(null);renderMatrix();$('download').disabled=true;
  $('energy-value').textContent=mathText('...');badge('RUNNING');
  $('status').textContent='Computing molecular integrals...';$('summary').textContent='Running locally in your browser.';
  document.querySelector('.summary-details').open=false;
  $('result-context').textContent=`${request.mode==='integrals'?'INTEGRALS':request.options.method.toUpperCase()} / ${basisLabel(request.basis)}`;
  $('job-progress').hidden=false;$('progress-fill').style.width='0%';
  $('progress-text').textContent=request.mode==='scan'?`0 / ${spec.count} scan points`:'Starting SCF';
  if(request.mode==='scan') {
    scan={spec,source:structuredClone(target),moleculeName:$('workspace-title').textContent,formula:$('formula').textContent,options:{...request.options},points:[],cancelled:false};
    selectedPoint=-1;scanRevision=draftRevision;renderScan();setTab('pes');
  } else setTab(TASKS.find(task=>task.id===request.mode)?.tab??'scf');
  setBusy(true);
  if(matchMedia('(max-width:850px)').matches)toggleSetup(false);
  // All validation, labels and computation use this one immutable snapshot.
  worker.postMessage(request);
}
function selectPoint(index,{restoreSettings=true}={}) {
  // An old plot remains inspectable after a scan, but must never change a
  // different running calculation's setup.
  if(busy&&jobMode!=='scan')return;
  const point=scan?.points[index];if(!point)return;
  const preserveDraft=scanRevision!==draftRevision;
  if(!busy){clearTimeout(previewTimer);++draftRevision;++jobId;activeRequest=null;}
  selectedPoint=index;
  $('geometry').value=toXYZ(point.molecule,'angstrom',`Scan point ${index+1}`);$('unit').value='angstrom';if(!preserveDraft)$('basis').value=scan.spec.basis;
  $('preset').value='custom';
  if(restoreSettings&&!busy&&!preserveDraft) {
    $('method').value=scan.options.method??'rhf';$('multiplicity').value=scan.options.multiplicity??1;$('guess').value=scan.options.guess??'core';
    $('charge').value=scan.options.charge;$('maxiter').value=scan.options.maxIterations;$('diis').checked=scan.options.diis;
    methodExplicit=true;
  }
  preview();
  if(point.result){
    showResult(point.result,{pointIndex:index});
    if(preserveDraft)$('status').textContent+=` Saved scan result; current setup ${$('method').value.toUpperCase()}/${basisLabel($('basis').value)} is preserved. Run to recalculate.`;
  }
  else {
    current=null;history=[];startHistory();drawSCF();renderOrbitals(null);renderMatrix();
    $('ao-order').textContent='';$('eri-selectors').replaceChildren();renderERI();
    $('results-title').textContent=`Scan point ${index+1} / ${scan.spec.count}`;
    $('metric-cycles').textContent='...';document.querySelector('.summary-details').open=true;
    badge('POINT FAILED','error');$('energy-value').textContent=mathText('...');$('status').textContent=`Scan point ${index+1} failed.`;$('summary').textContent=point.error;$('download').disabled=true;
  }
  renderScan();previewPoint(-1);
  if(!busy){if(!preserveDraft)scanRevision=draftRevision;persistSettings();}
}
function previewPoint(index) {
  if(busy&&jobMode!=='scan')return;
  const point=scan?.points[index];
  if(point){viewer.setMolecule(point.molecule);$('pes-preview').textContent=mathText(`Preview: point ${index+1} · ${point.coordinate.toFixed(4)} ${scan.spec.unit==='degree'?'°':scan.spec.unit==='angstrom'?'Å':'bohr'} · ${point.status}`);}
  else {if(molecule)viewer.setMolecule(molecule);$('pes-preview').textContent='Hover a point to preview its geometry. Click to inspect its result.';}
}
function renderScan() {
  const points=scan?.points??[];$('scan-count-badge').textContent=points.length;
  for(const id of ['scan-csv','scan-xyz','scan-json'])$(id).disabled=!points.length;
  document.querySelectorAll('[data-pes-export]').forEach(button=>button.disabled=!points.some(p=>Number.isFinite(p.energy)));
  if(!scan){plot($('pes-chart'),[],{empty:'Your potential energy curve will appear here.'});$('scan-point-table').replaceChildren();return;}
  const figure=points.some(p=>Number.isFinite(p.energy))?pesFigureModel(scan,{unit:$('pes-energy-unit').value,relative:$('pes-relative').checked}):null;
  if(figure)plot($('pes-chart'),figure.points,{xLabel:figure.xLabel,yLabel:figure.yLabel,onSelect:selectPoint,onPreview:previewPoint,selected:selectedPoint});
  else plot($('pes-chart'),[],{empty:points.length?'No numerical energies are available yet.':'Scan points will appear as calculations complete.'});
  $('pes-note').textContent=`${basisLabel(scan.spec.basis)} · unrelaxed ${(scan.spec.method??'rhf').toUpperCase()} · ${scan.spec.coordinate??'bond'} · ${points.length}/${scan.spec.count} points${scan.cancelled?' · stopped':''}. ${figure?.relative?'Zero at the lowest '+(points.some(p=>p.converged)?'converged':'available, unconverged')+' sampled point.':''} Failed points are not connected.`;
  const meta=COORDINATES[scan.spec.coordinate??'bond'];
  const body=table($('scan-point-table'),['Point',`${meta.symbol} / ${scan.spec.unit==='degree'?'°':scan.spec.unit==='angstrom'?'Å':'bohr'}`,'Energy / Eₕ','Status','Cycles'],points.map(p=>[p.index+1,p.coordinate.toFixed(6),p.energy===null?'Not available':p.energy.toFixed(12),p.status,p.iterations]));
  [...body.rows].forEach((row,index)=>{row.className=`clickable ${selectedPoint===index?'selected':''}`;row.tabIndex=0;row.setAttribute('aria-label',`Inspect scan point ${index+1}`);row.addEventListener('click',()=>selectPoint(index));row.addEventListener('pointerenter',()=>previewPoint(index));row.addEventListener('pointerleave',()=>previewPoint(-1));row.addEventListener('keydown',e=>{if(e.key==='Enter')selectPoint(index);});});
  if(scan?.summary&&scan.listId==null)scan.listId=rememberCompletedJob(current,{kind:'scan'});
}
worker.onmessage=({data})=>{
  if(data.type==='ready') {
    if(data.protocol!==WORKER_PROTOCOL||data.buildId!==BUILD_ID) {
      ready=false;current=null;worker.terminate();setBusy(false);$('energy-value').textContent=mathText('...');badge('ENGINE VERSION ERROR','error');
      $('status').textContent='Page and worker versions differ. Reload the page.';
      $('summary').textContent='Restart the local server and reload the page to load one consistent build.';
      document.querySelector('.summary-details').open=true;return;
    }
    ready=true;$('engine-state').textContent='Local engine ready';$('engine-dot').classList.add('ready');
    preview();setBusy(false);if(draftRevision===0)run();return;
  }
  if(data.type==='fatal') {$('engine-state').textContent='Engine unavailable';ready=false;setBusy(false);badge('ENGINE ERROR','error');$('status').textContent='Engine unavailable.';$('summary').textContent=data.message;return;}
  if(data.id!==jobId||!activeRequest)return;
  if(data.protocol!==WORKER_PROTOCOL||data.buildId!==BUILD_ID) {
    activeRequest=null;current=null;ready=false;worker.terminate();setBusy(false);$('energy-value').textContent=mathText('...');
    badge('ENGINE VERSION ERROR','error');$('status').textContent='Page and worker versions differ. Reload the page.';
    $('summary').textContent='Calculation discarded because it came from a different engine build.';
    document.querySelector('.summary-details').open=true;return;
  }
  if(data.type==='post-scf') {
    const p=data.progress;$('status').textContent=`MP2 correlation: occupied block ${p.completed} / ${p.occupied}`;
    $('progress-text').textContent=`${data.index===undefined?'':`Point ${data.index+1} · `}MP2 ${p.completed} / ${p.occupied} occupied orbitals`;return;
  }
  if(data.type==='iteration') {
    $('energy-label').textContent=activeRequest.options.method==='mp2'?'RHF REFERENCE IN PROGRESS':'SCF ENERGY IN PROGRESS';
    if(jobMode==='scan') {
      $('progress-text').textContent=`Point ${data.index+1} / ${scan.spec.count} · SCF cycle ${data.entry.iteration}`;
      if(!current){$('status').textContent=`SCF iteration ${data.entry.iteration} at scan point ${data.index+1}`;$('energy-value').textContent=mathText(data.entry.energy.toFixed(10));}
    }else{
      history.push(data.entry);appendIteration(data.entry);drawSCF();
      $('status').textContent=mathText(`SCF iteration ${data.entry.iteration}: E = ${data.entry.energy.toFixed(10)} Eₕ`);$('energy-value').textContent=mathText(data.entry.energy.toFixed(10));
      $('progress-text').textContent=`SCF cycle ${data.entry.iteration}`;
    }return;
  }
  if(data.type==='point') {
    const follow=selectedPoint<0||selectedPoint===scan.points.length-1;
    if(data.point.result) {
      try {assertResultMatches(activeRequest,data.point.result);}
      catch(error) {worker.postMessage({type:'cancel',id:jobId});activeRequest=null;setBusy(false);badge('ENGINE ERROR','error');$('status').textContent=error.message;$('summary').textContent=error.message;return;}
    }
    scan.points.push(data.point);scan.spec=data.spec;
    $('progress-fill').style.width=`${100*scan.points.length/scan.spec.count}%`;
    if(follow)selectPoint(scan.points.length-1,{restoreSettings:false});else renderScan();return;
  }
  if(data.type==='scan-complete') {
    scan.cancelled=data.summary.cancelled;scan.summary=data.summary;activeRequest=null;setBusy(false);renderScan();persistSettings();
    $('progress-text').textContent=`${scan.cancelled?'Stopped':'Complete'} · ${scan.points.length} / ${scan.spec.count} points`;
    if(!scan.points.length){badge('CANCELLED','warning');$('status').textContent='Scan cancelled before a point completed.';}
    return;
  }
  setBusy(false);$('job-progress').hidden=true;
  if(data.type==='result'){
    const request=activeRequest;activeRequest=null;
    try {assertResultMatches(request,data.result);showResult(data.result);}
    catch(error) {current=null;$('download').disabled=true;badge('ENGINE ERROR','error');$('energy-value').textContent=mathText('...');$('status').textContent=error.message;$('summary').textContent=error.message;document.querySelector('.summary-details').open=true;}
    return;
  }
  activeRequest=null;
  current=null;$('download').disabled=true;$('energy-value').textContent=mathText('...');$('metric-cycles').textContent='...';
  document.querySelector('.summary-details').open=data.type!=='cancelled';
  badge(data.type==='cancelled'?'CANCELLED':'INPUT ERROR',data.type==='cancelled'?'warning':'error');
  $('status').textContent=data.type==='cancelled'?'Calculation cancelled.':'Input or calculation error.';$('summary').textContent=data.message??'You can adjust the settings and run again.';
};
worker.onerror=event=>{ready=false;setBusy(false);badge('WORKER ERROR','error');$('engine-state').textContent='Engine unavailable';$('status').textContent='The calculation worker failed.';$('summary').textContent=`${event.message}\nReload the page to restart.`;};
function save(name,content,type='text/plain') {
  const url=URL.createObjectURL(new Blob([content],{type})),a=document.createElement('a');a.href=url;a.download=name;a.click();setTimeout(()=>URL.revokeObjectURL(url),1000);
}
$('controls').addEventListener('submit',run);
$('controls').addEventListener('invalid',event=>{
  const details=event.target.closest('details');
  if(details)details.open=true;
},true);
$('cancel').addEventListener('click',()=>{worker.postMessage({type:'cancel',id:jobId});$('cancel').disabled=true;$('cancel').textContent='Stopping...';$('progress-text').textContent='Stopping after the current integral or SCF step...';});
function settingsChanged({reset=false}={}) {
  invalidate();preview({reset});persistSettings();
}
$('controls').addEventListener('input',event=>{
  if(!SETTINGS_FIELDS.includes(event.target.id))return;
  if(busy||['SELECT','BUTTON'].includes(event.target.tagName)||event.target.id==='xyz-file')return;
  if(event.target.id==='geometry')$('preset').value='custom';
  invalidate();persistSettings();
  if(event.target.id.startsWith('scan-'))selection();
});
$('controls').addEventListener('change',event=>{
  if(!SETTINGS_FIELDS.includes(event.target.id))return;
  if(busy||event.target.id==='xyz-file')return;
  const id=event.target.id;
  if(id==='method')methodExplicit=true;
  if(id==='preset') {
    const preset=$('preset').value,settings=EXAMPLE_SETTINGS[preset];
    if(!EXAMPLES[preset])return;
    $('geometry').value=EXAMPLES[preset];$('unit').value='angstrom';
    $('charge').value=settings?.charge??0;$('multiplicity').value=settings?.multiplicity??1;
    $('method').value=presetMethod($('method').value,settings?.method??'rhf',methodExplicit);
    if(settings?.method==='uhf')$('guess').value='core';
    $('scan-coordinate').value='bond';$('scan-unit').value='angstrom';$('scan-a').value='1';$('scan-b').value='2';$('scan-c').value='';$('scan-d').value='';$('scan-moving').value='2';
  }
  settingsChanged({reset:id==='preset'});
  if(id==='mode')$('task-popover').open=false;
});
$('water').addEventListener('click',()=>{
  if(busy)return;
  methodExplicit=false;$('controls').reset();$('scan-coordinate').value='bond';$('scan-unit').value='angstrom';$('geometry').value=WATER_XYZ;viewer.setLabels($('atom-labels').checked);
  settingsChanged({reset:true});run();
});
$('xyz-file').addEventListener('change',async event=>{
  const file=event.target.files[0];if(!file)return;
  if(file.size>1000000){$('input-message').textContent='Use an XYZ file smaller than 1 MB.';return;}
  const revision=draftRevision,text=await file.text();if(busy||revision!==draftRevision)return;
  $('geometry').value=text;$('preset').value='custom';$('unit').value='angstrom';settingsChanged({reset:true});event.target.value='';
});
$('atom-labels').addEventListener('change',()=>viewer.setLabels($('atom-labels').checked));
$('fit-view').addEventListener('click',()=>viewer.reset());
$('viewer-png').addEventListener('click',()=>{const uri=viewer.png();if(!uri){$('input-message').textContent='Image export needs a working WebGL viewer.';return;}const a=document.createElement('a');a.href=uri;a.download='molecule.png';a.click();});
$('geometry-download').addEventListener('click',()=>{if(molecule)save('geometry.xyz',toXYZ(molecule,'angstrom',`${basisLabel(molecule.basis)} / exported geometry`));});
$('matrix-kind').addEventListener('change',renderMatrix);
$('pes-energy-unit').addEventListener('change',renderScan);
$('pes-relative').addEventListener('change',renderScan);
document.querySelectorAll('[data-pes-export]').forEach(button=>button.addEventListener('click',async()=>{
  if(!scan?.points.some(p=>Number.isFinite(p.energy)))return;
  button.disabled=true;
  try{await exportPESFigure(pesFigureModel(scan,{unit:$('pes-energy-unit').value,relative:$('pes-relative').checked}),button.dataset.pesExport);}
  catch(error){$('pes-preview').textContent=`Figure export failed: ${error.message}`;}
  finally{renderScan();}
}));
$('download').addEventListener('click',()=>{if(current)save('calculation.json',JSON.stringify(serializable(current),null,2),'application/json');});
$('scan-csv').addEventListener('click',()=>save('pes-scan.csv',scanCSV(scan),'text/csv'));
$('scan-xyz').addEventListener('click',()=>save('pes-scan.xyz',scanXYZ(scan)));
$('scan-json').addEventListener('click',()=>save('pes-scan.json',JSON.stringify({...scan,points:scan.points.map(p=>({...p,...(p.result?{result:serializable(p.result)}:{})}))},null,2),'application/json'));
document.addEventListener('keydown',event=>{
  if((event.metaKey||event.ctrlKey)&&event.key==='Enter'&&!$('workspace-guide').open&&!busy){event.preventDefault();run();}
});
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
globalThis.addEventListener('pageshow',event=>{if(event.persisted&&!busy){invalidate();preview();}});

function rememberCompletedJob(result,{kind}={}) {
  const list=$('job-list');
  if(!list) return 0;
  const scanJob=kind==='scan';
  if(!result&&!scanJob) return 0;
  const r=result?scfResult(result):null;
  const method=scanJob?(scan?.options?.method??result?.method??'rhf').toUpperCase():(result.method??r?.method??'INTEGRALS');
  const basis=basisLabel(scanJob?scan?.spec?.basis??result?.molecule?.basis??$('basis').value:result.molecule.basis);
  const title=$('workspace-title').textContent||'Custom molecule';
  const status=scanJob?(scan.cancelled?'STOPPED':`${scan.points.length} pts`):($('result-badge').textContent||'DONE');
  const energy=scanJob?`${scan.points.filter(p=>p.converged).length} converged`:($('energy-value').textContent||'');
  const id=++listedJobSeq;
  listedJobs.unshift({
    id,kind:scanJob?'scan':result.checks?'integrals':'energy',title,method,basis,status,energy,
    fields:readSettings($('controls')),methodExplicit,result:result??null,scan:scanJob?scan:null,
    selectedPoint:scanJob?selectedPoint:-1
  });
  if(listedJobs.length>24) listedJobs.length=24;
  activeListedJob=id;
  renderJobList();
  return id;
}
function renderJobList() {
  const list=$('job-list');
  if(!list) return;
  $('job-count').textContent=String(listedJobs.length);
  list.replaceChildren();
  for(const job of listedJobs) {
    const item=document.createElement('button');
    item.type='button';
    item.className=`job-item${job.id===activeListedJob?' active':''}`;
    item.innerHTML=`<span class="job-status">${job.status}</span><span class="job-name"></span><span class="job-meta"></span>`;
    item.querySelector('.job-name').textContent=job.title;
    item.querySelector('.job-meta').textContent=`${job.method}/${job.basis} · ${job.energy}`;
    item.addEventListener('click',()=>openListedJob(job.id));
    const row=document.createElement('li');
    row.append(item);
    list.append(row);
  }
}
function openListedJob(id) {
  if(busy) return;
  const job=listedJobs.find(entry=>entry.id===id);
  if(!job) return;
  activeListedJob=id;
  if(job.fields) {
    for(const [key,value] of Object.entries(job.fields)) {
      if(key==='diis') $('diis').checked=value;
      else if($(key)&&value!=null) $(key).value=value;
    }
    methodExplicit=!!job.methodExplicit;
  }
  preview();
  if(job.kind==='scan'&&job.scan) {
    scan=job.scan;selectedPoint=job.selectedPoint??scan.points.length-1;renderScan();
    const point=scan.points[selectedPoint];
    if(point?.result) showResult(point.result,{pointIndex:selectedPoint,skipRemember:true});
  } else if(job.result) showResult(job.result,{skipRemember:true});
  persistSettings();
  renderJobList();
}
renderJobList();
