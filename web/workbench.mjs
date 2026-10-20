// SPDX-License-Identifier: AGPL-3.0-only
// New tasks, methods and tabs register here without changing the workspace grid.
export const TASKS=[
  {id:'rhf',label:'Single-point energy',tab:'scf',solver:true},
  {id:'scan',label:'Bond-length scan',tab:'pes',solver:true},
  {id:'integrals',label:'Molecular integrals',tab:'matrices',solver:false}
];
export const TABS=[
  {id:'scf',label:'SCF convergence'}, {id:'orbitals',label:'Orbitals'},
  {id:'matrices',label:'AO matrices'}, {id:'pes',label:'Potential energy',count:true},
  {id:'data',label:'Geometry'}
];
export const METHODS={rhf:{groups:['scf']},uhf:{groups:['scf','uhf']},mp2:{groups:['scf','mp2']}};
const $=id=>document.getElementById(id);
export function closeMenus(except=null) {
  document.querySelectorAll('.popover[open]').forEach(panel=>{if(panel!==except)panel.open=false;});
}
export function toggleSetup(open) {
  if(!open)closeMenus();
  $('coordinates-editor').open=!!open;
}
export function focusViewer(enabled) {
  document.body.classList.toggle('focus-view',enabled);
  $('focus-view').setAttribute('aria-pressed',String(enabled));
  $('focus-view').title=enabled?'Exit expanded viewer':'Expand molecule viewer';
}
export function initializeWorkbench() {
  $('mode').replaceChildren(...TASKS.map(task=>new Option(task.label,task.id)));
  const tabs=document.querySelector('.result-tabs');
  for(const [index,tab] of TABS.entries()) {
    const button=document.createElement('button');
    button.type='button';button.id=`tab-${tab.id}`;button.dataset.tab=tab.id;
    button.setAttribute('role','tab');button.setAttribute('aria-controls',`panel-${tab.id}`);
    button.setAttribute('aria-selected',String(index===0));button.tabIndex=index===0?0:-1;
    button.textContent=tab.label;
    if(tab.count){const count=document.createElement('span');count.id='scan-count-badge';count.textContent='0';button.append(count);}
    tabs.append(button);
    $(`panel-${tab.id}`).setAttribute('aria-labelledby',button.id);
  }
  function syncRecipe() {
    const task=TASKS.find(t=>t.id===$('mode').value);
    $('task-label').textContent=task.label;
    const mult=$('method').value==='uhf'?Number($('multiplicity').value):1;
    const spin={1:'Singlet',2:'Doublet',3:'Triplet',4:'Quartet',5:'Quintet'}[mult]??`2S+1: ${mult}`;
    $('spin-label').textContent=`Q: ${$('charge').value} · ${task.solver?spin:'Integrals'}`;
    document.querySelectorAll('[data-solver-group]').forEach(group=>group.hidden=!task.solver||!METHODS[$('method').value].groups.includes(group.dataset.solverGroup));
    $('integrals-solver-note').hidden=task.solver;
  }
  $('controls').addEventListener('input',syncRecipe);
  $('controls').addEventListener('change',syncRecipe);
  // Result inspection and restored jobs also change fields programmatically.
  new MutationObserver(syncRecipe).observe($('molecule-state'),{childList:true});
  document.querySelectorAll('.popover').forEach(panel=>panel.addEventListener('toggle',()=>{
    panel.querySelector('summary').setAttribute('aria-expanded',String(panel.open));
    if(panel.open)closeMenus(panel);
  }));
  document.addEventListener('click',event=>{
    if(!event.target.closest('.popover'))closeMenus();
    if(!$('jobs-drawer').hidden&&!event.target.closest('#jobs-drawer, #jobs-toggle'))jobs(false);
  });
  function jobs(open) {
    $('jobs-drawer').hidden=!open;$('jobs-toggle').setAttribute('aria-expanded',String(open));
    if(open){closeMenus();$('jobs-close').focus();}
  }
  $('jobs-toggle').addEventListener('click',()=>jobs($('jobs-drawer').hidden));
  $('jobs-close').addEventListener('click',()=>{jobs(false);$('jobs-toggle').focus();});
  $('xyz-toggle').addEventListener('click',()=>{
    $('coordinates-editor').open=!$('coordinates-editor').open;
    if($('coordinates-editor').open)$('geometry').focus();
  });
  $('coordinates-editor').addEventListener('toggle',()=>{
    $('xyz-toggle').setAttribute('aria-expanded',String($('coordinates-editor').open));
    if(!$('coordinates-editor').open&&$('coordinates-editor').contains(document.activeElement))$('xyz-toggle').focus();
  });
  $('focus-view').addEventListener('click',()=>focusViewer(!document.body.classList.contains('focus-view')));
  $('help-open').addEventListener('click',()=>{closeMenus();$('workspace-guide').showModal();});
  for(const id of ['help-close','help-done'])$(id).addEventListener('click',()=>$('workspace-guide').close());
  $('geometry-copy').addEventListener('click',async()=>{
    try{await navigator.clipboard.writeText($('geometry-output').textContent);$('geometry-copy').textContent='Copied';}
    catch{$('geometry-copy').textContent='Copy unavailable';}
    setTimeout(()=>$('geometry-copy').textContent='Copy XYZ',1800);
  });
  document.addEventListener('keydown',event=>{
    if(event.key!=='Escape'||$('workspace-guide').open)return;
    const open=document.querySelector('.popover[open]');
    if(open){open.open=false;open.querySelector('summary').focus();return;}
    if(!$('jobs-drawer').hidden){jobs(false);$('jobs-toggle').focus();return;}
    if($('coordinates-editor').open){$('coordinates-editor').open=false;$('xyz-toggle').focus();return;}
    focusViewer(false);
  });
  const handle=$('split-handle'),workbench=$('workbench');let dragging=false;
  function split(value) {
    const clamped=Math.max(35,Math.min(70,value));
    workbench.style.setProperty('--stage-width',`${clamped}%`);handle.setAttribute('aria-valuenow',String(Math.round(clamped)));
  }
  handle.addEventListener('pointerdown',event=>{
    if(event.button!==0)return;
    event.preventDefault();dragging=true;handle.setPointerCapture(event.pointerId);workbench.classList.add('resizing');
  });
  handle.addEventListener('pointermove',event=>{
    if(!dragging)return;
    const box=workbench.getBoundingClientRect();split(100*(event.clientX-box.left)/box.width);
  });
  function stop(){dragging=false;workbench.classList.remove('resizing');}
  handle.addEventListener('pointerup',stop);handle.addEventListener('pointercancel',stop);handle.addEventListener('lostpointercapture',stop);
  handle.addEventListener('keydown',event=>{
    const v=Number(handle.getAttribute('aria-valuenow'));
    const next={ArrowLeft:v-2,ArrowRight:v+2,Home:35,End:70}[event.key];
    if(next!==undefined){event.preventDefault();split(next);}
  });
  handle.addEventListener('dblclick',()=>split(55));
  document.addEventListener('lorentzian:orbital-select',({detail})=>{
    const slot=$('orbital-overlay-slot');slot.hidden=!detail;slot.replaceChildren();if(!detail)return;
    const label=document.createElement('span');label.textContent=`${detail.spin?detail.spin+' ':''}MO ${detail.index+1} · ${detail.energy.toFixed(6)} Eₕ`;
    const note=document.createElement('small');note.textContent='Orbital selection';slot.append(label,note);
  });
  syncRecipe();
}
