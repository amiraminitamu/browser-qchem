// SPDX-License-Identifier: AGPL-3.0-only
// Requires jsdom as a development-only dependency. No layout or WebGL assertions.
import {createRequire} from 'node:module';
const require=createRequire(process.env.LORENTZIAN_QA_PACKAGE??import.meta.url);
const {JSDOM}=require('jsdom');
import fs from 'node:fs';
import assert from 'node:assert/strict';
import {Worker as NodeWorker} from 'node:worker_threads';
import {pathToFileURL,fileURLToPath} from 'node:url';
const root=fileURLToPath(new URL('../',import.meta.url));
const dom=new JSDOM(fs.readFileSync(root+'/build/wasm/dist/index.html','utf8'),{url:'https://lorentzian.test/',pretendToBeVisual:true});
const {window}=dom;
for(const key of ['document','MutationObserver','CustomEvent','Option','HTMLElement','HTMLCanvasElement','Event'])globalThis[key]=window[key];
globalThis.addEventListener=window.addEventListener.bind(window);
globalThis.matchMedia=()=>({matches:false});
globalThis.ResizeObserver=class{observe(){}};
const workers=[],errors=[];
// Adapt browser URL loading to Node filesystem URLs; retain the actual worker
// dispatch, protocol handshake, solver code, cancellation and streamed messages.
async function bootstrap(){
  const {parentPort,workerData}=await import('node:worker_threads');
  const fs=await import('node:fs');
  globalThis.self={postMessage:data=>parentPort.postMessage(data)};
  parentPort.on('message',data=>self.onmessage({data}));
  let source=fs.readFileSync(new URL(workerData.url),'utf8');
  source=source.replace("new URL(path,prefix||import.meta.url)","new URL(path,prefix?.startsWith('/')?'file://'+prefix:prefix||import.meta.url)").replace("url.searchParams.set('v',BUILD_ID);","if(url.protocol!=='file:')url.searchParams.set('v',BUILD_ID);");
  source=source.replace(/from '(\.[^']+)'/g,(_match,path)=>"from '"+new URL(path,workerData.url).href+"'");
  await import('data:text/javascript,'+encodeURIComponent(source));
}
class WorkerShim {
  constructor(url){this.worker=new NodeWorker(new URL('data:text/javascript,'+encodeURIComponent('('+bootstrap.toString()+')()')),{workerData:{url:url.href}});workers.push(this.worker);
    this.worker.on('message',data=>this.onmessage?.({data}));this.worker.on('error',error=>{errors.push(error);this.onerror?.(error);});}
  postMessage(data){this.worker.postMessage(data);}
  terminate(){this.worker.terminate();}
}
globalThis.Worker=WorkerShim;
const $=id=>document.getElementById(id);
const event=(id,type='change')=>$(id).dispatchEvent(new window.Event(type,{bubbles:true}));
const set=(id,value)=>{$(id).value=value;event(id);};
const waitFor=async(fn,label)=>{const end=Date.now()+15000;while(!fn()){if(Date.now()>end)throw new Error(label+': '+$('status').textContent+' / '+$('input-message').textContent+' / '+$('summary').textContent);await new Promise(r=>setTimeout(r,15));}};
try{
  await import(pathToFileURL(root+'/build/wasm/dist/app.mjs'));
  await waitFor(()=>$('status').textContent.startsWith('CONVERGED'),'initial water');
  assert.equal($('scf-integrals').value,'df');assert.match($('method-details').textContent,/Density fitted \(DF-JK\)/);
  await waitFor(()=>$('eri-value').textContent.includes('exact, on demand'),'DF exact ERI inspector');
  const dfERI=$('eri-value').textContent.split(' Eₕ')[0];
  set('scf-integrals','exact');$('calculate').click();await waitFor(()=>$('status').textContent.startsWith('CONVERGED'),'exact water');
  assert.match($('method-details').textContent,/Exact four-center/);assert.equal($('eri-value').textContent,dfERI+' Eₕ');
  set('scf-integrals','df');$('calculate').click();await waitFor(()=>$('status').textContent.startsWith('CONVERGED'),'DF water again');
  assert.equal($('energy-value').textContent[0],'−');assert.equal($('scan-a').value,'1');assert.equal($('scan-b').value,'2');
  assert.equal($('scan-controls').hidden,true);assert.equal($('selection-info').hidden,true);
  set('mode','scan');document.querySelector('[data-coordinate=angle]').click();
  assert.equal($('scan-unit').value,'degree');
  for(const [id,value] of [['a','2'],['b','1'],['c','3']])set('scan-'+id,value);
  assert.equal($('scan-moving').value,'3');assert.match($('scan-current').textContent,/104\./);
  set('scan-start','90');set('scan-end','120');set('scan-count','3');
  $('calculate').click();
  await waitFor(()=>$('progress-text').textContent.startsWith('Complete'),'angle scan');
  assert.equal($('scan-point-table').querySelectorAll('tbody tr').length,3);
  const row=$('scan-point-table').querySelectorAll('tbody tr')[0];row.dispatchEvent(new window.Event('pointerenter'));
  assert.match($('pes-preview').textContent,/90\.0000 °/);row.dispatchEvent(new window.Event('pointerleave'));
  row.click();assert.equal($('preset').value,'custom');
  set('pes-energy-unit','ev');assert.match($('pes-chart').textContent,/eV/);
  set('mode','rhf');assert.equal($('scan-controls').hidden,true);assert.equal($('selection-info').hidden,true);
  set('method','mp2');$('calculate').click();await waitFor(()=>$('status').textContent.startsWith('CONVERGED: MP2'),'MP2 after scan');
  assert.match($('summary').textContent,/Total MP2 energy/);
  set('mode','scan');document.querySelector('[data-coordinate=dihedral]').click();
  $('geometry').value='H -.3 .9 0\nO 0 0 0\nO 1.45 0 0\nH 1.75 0 .9';event('geometry');
  for(const [id,value] of [['a','1'],['b','2'],['c','3'],['d','4']])set('scan-'+id,value);
  assert.equal($('scan-moving').value,'3, 4');set('scan-start','30');set('scan-end','150');set('scan-count','3');set('method','rhf');$('calculate').click();
  await waitFor(()=>$('progress-text').textContent.startsWith('Complete'),'dihedral scan');
  assert.equal($('scan-point-table').querySelectorAll('tbody tr').length,3);assert.match($('pes-chart').textContent,/φ/);
  $('scan-clear').click();assert.equal($('scan-a').value,'');assert.equal($('scan-d').value,'');
  $('water').click();await waitFor(()=>$('status').textContent.startsWith('CONVERGED'),'full reset');
  assert.equal($('scan-coordinate').value,'bond');assert.equal($('scan-a').value,'1');assert.equal($('scan-unit').value,'angstrom');
  set('method','df-rmp2');$('calculate').click();await waitFor(()=>$('status').textContent.startsWith('CONVERGED: DF-RMP2'),'DF-RMP2 water');
  assert.equal($('uhf-controls').hidden,true);assert.match($('method-details').textContent,/cc-pVDZ-RI/);assert.match($('method-details').textContent,/MP2 correlation/);
  assert.equal(document.querySelector('#panel-scf h3').textContent,'RHF reference convergence');
  set('method','df-ump2');set('preset','h2plus');$('calculate').click();await waitFor(()=>$('status').textContent.startsWith('CONVERGED: DF-UMP2'),'DF-UMP2 ion');
  assert.equal($('method').value,'df-ump2');assert.equal($('uhf-controls').hidden,false);assert.match($('method-details').textContent,/Spin contamination/);assert.match($('method-details').textContent,/MP2 correlation/);
  assert.equal(document.querySelector('#panel-scf h3').textContent,'UHF reference convergence');
  set('mode','scan');set('scan-start','.8');set('scan-end','1.2');set('scan-count','3');$('calculate').click();
  await waitFor(()=>$('progress-text').textContent.startsWith('Complete'),'DF-UMP2 ion scan');assert.match($('pes-note').textContent,/DF-UMP2/);
  assert.equal($('scan-point-table').querySelectorAll('tbody tr').length,3);
  assert.ok(!errors.length,errors.map(String).join('\n'));
  console.log('DOM + actual worker/WASM: initial RHF, mode gating, angle and torsion UI scans, fragments, plot preview, unit conversion, MP2, both DF methods, spin controls and DF PES PASS');
}finally{for(const worker of workers)await worker.terminate();window.close();}
