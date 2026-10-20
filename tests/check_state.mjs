// SPDX-License-Identifier: AGPL-3.0-only
import fs from 'node:fs';
import {mathText} from '../web/format.mjs';
import {TASKS} from '../web/workbench.mjs';

function check(value,message) {if(!value)throw new Error(message);}
function equal(a,b,message) {check(Object.is(a,b),message+': '+JSON.stringify(a)+' != '+JSON.stringify(b));}
function exported(source,names) {return Function(source.replace(/^export /gm,'')+'\nreturn {'+names.join(',')+'};')();}
function functionBlock(source,start,end) {const a=source.indexOf(start),b=source.indexOf(end,a);check(a>=0&&b>a,'Controller test source anchors');return source.slice(a,b);}

export function runStateTests(appSource,settingsSource) {
  const settings=exported(settingsSource,['readSettings','restoreSettings','saveSettings','presetMethod','captureRequest','assertResultMatches','SETTINGS_FIELDS']);
  const defaults={preset:'water',geometry:'O 0 0 0\nH 0 0 1\nH 0 1 0',unit:'angstrom',basis:'sto-3g',method:'rhf',mode:'rhf',charge:'0',multiplicity:'1',guess:'mix',maxiter:'100',diis:true,'scan-coordinate':'bond','scan-c':'','scan-d':'','scan-a':'1','scan-b':'2','scan-moving':'2','scan-start':'.7','scan-end':'1.4','scan-count':'3','scan-unit':'angstrom'};
  function controller(saved) {
    const nodes=new Map(),queue=new Map(),memory=new Map(saved??[]);let clock=0;
    function node(id) {
      if(nodes.has(id))return nodes.get(id);
      const listeners={},classes=new Set();
      const n={id,value:String(defaults[id]??''),checked:id==='diis',disabled:false,tagName:['preset','unit','basis','method','mode','guess','scan-unit'].includes(id)?'SELECT':'INPUT',textContent:'',hidden:false,open:false,style:{},
        classList:{add:v=>classes.add(v),contains:v=>classes.has(v),toggle:(v,on)=>on?classes.add(v):classes.delete(v)},
        addEventListener:(type,fn)=>(listeners[type]??=[]).push(fn),
        emit(type){const event={target:n,preventDefault(){}};for(const fn of listeners[type]??[])fn(event);if(id!=='controls')for(const fn of node('controls').listeners[type]??[])fn(event);},
        listeners,replaceChildren(){},querySelector:selector=>node(selector.replace(/^#/,'')),querySelectorAll:()=>[],
        closest:selector=>selector==='#scan-controls'&&id.startsWith('scan-')?node('scan-controls'):null,
        checkValidity:()=>true,reportValidity:()=>true,scrollIntoView(){},focus(){},setAttribute(){},contains:()=>false};
      if(id==='controls') {
        n.querySelectorAll=()=>settings.SETTINGS_FIELDS.map(node);
        n.reset=()=>{for(const [key,value]of Object.entries(defaults)){if(key==='diis')node(key).checked=value;else node(key).value=value;}};
      }
      let nativeValue=n.value;Object.defineProperty(n,'value',{get:()=>nativeValue,set:value=>{nativeValue=String(value);}});
      nodes.set(id,n);return n;
    }
    const form=node('controls'),storage={getItem:k=>memory.get(k)??null,setItem:(k,v)=>memory.set(k,v)};
    const explicit=settings.restoreSettings(form,storage);
    const sent=[],shown=[];
    const env={$:node,document:{querySelector:selector=>node(selector),querySelectorAll:()=>[],dispatchEvent(){}},worker:{postMessage:r=>sent.push(r),terminate(){}},
      BUILD_ID:'test-build',WORKER_PROTOCOL:2,EXAMPLES:{water:defaults.geometry,hehplus:'He 0 0 0\nH 0 0 .8',h2plus:'H 0 0 0\nH 0 0 1.06'},EXAMPLE_SETTINGS:{hehplus:{method:'rhf',charge:1,multiplicity:1},h2plus:{method:'uhf',charge:1,multiplicity:2}},WATER_XYZ:defaults.geometry,
      ...settings,TASKS,mathText,viewer:{setLabels(){},setMolecule(){}},structuredClone:v=>JSON.parse(JSON.stringify(v)),Number,Object,Array,Math,Set,
      setTimeout:fn=>{queue.set(++clock,fn);return clock;},clearTimeout:id=>queue.delete(id),
      parseGeometry:(geometry,unit,basis)=>({geometry,unit,basis,atoms:[{charge:1,xyz:[0,0,0]},{charge:1,xyz:[0,0,2]}]}),
      validateElectronicOptions:(_m,o)=>{if(!['rhf','uhf','mp2'].includes(o.method))throw Error('Invalid method');return o;},
      scanSpec:(_m,c)=>({...c,basis:node('basis').value}),basisLabel:b=>b,
      matchMedia:()=>({matches:false}),referenceContext(){},plot(){},startHistory(){},drawSCF(){},renderOrbitals(){},renderMatrix(){},setTab(){},renderScan(){},selection(){},toggleSetup(){},toXYZ:()=> 'H 0 0 0\nH 0 0 2',
      showResult:r=>shown.push(r),badge:(text)=>{node('result-badge').textContent=text;},renderERI(){}};
    const setup=functionBlock(appSource,'function setBusy(value)','function badge(');
    const disabled=appSource.includes('function setDisabled(')?functionBlock(appSource,'function setDisabled(','function setBusy(value)'):'';
    const invalidate=functionBlock(appSource,'function invalidate()','function table(');
    const run=functionBlock(appSource,'function run(event)','function selectPoint(');
    const point=functionBlock(appSource,'function selectPoint(','function renderScan()');
    const events=functionBlock(appSource,"$('controls').addEventListener('input',","$('xyz-file').addEventListener");
    const changed=appSource.includes('function settingsChanged(')?functionBlock(appSource,'function settingsChanged(',"$('controls').addEventListener('input',"):'';
    const message=functionBlock(appSource,'worker.onmessage=','worker.onerror=');
    const options=appSource.includes('function options()')?functionBlock(appSource,'function options()','function run(event)'):'';
    const source='with(env){let busy=false,ready=true,valid=true,jobId=0,jobMode="rhf",molecule={atoms:[]},current=null,scan=null,selectedPoint=-1,history=[],previewTimer=null,activeRequest=null,draftRevision=0,scanRevision=0,methodExplicit=initialExplicit,pair=[0,1],selectionExplicit=false;'+
      'function persistSettings(){saveSettings($("controls"),storage,methodExplicit);}'+
      'function preview(){clearTimeout(previewTimer);if(typeof syncControls==="function")syncControls();else if(!busy)setBusy(false);}'+
      disabled+setup+invalidate+options+run+point+changed+events+message+
      ';return {run,selectPoint,setBusy,preview,flush(){const pending=[...queue.values()];queue.clear();pending.forEach(f=>f());},setScan(value){scan=value;scanRevision=draftRevision;},snapshot(){return {busy,jobId,activeRequest,methodExplicit};}};}';
    const api=Function('env','storage','initialExplicit','queue',source)(env,storage,explicit,queue);
    return {api,node,form,storage,memory,sent,shown,queue,env};
  }
  const reports=[];
  const c=controller();
  c.node('method').value='mp2';c.node('method').emit('change');
  c.node('preset').value='hehplus';c.node('preset').emit('change');
  equal(c.node('method').value,'mp2','Explicit MP2 must survive a closed-shell ion preset');
  c.api.run();equal(c.sent.at(-1).options.method,'mp2','Submitted method after preset');
  check(Object.isFrozen(c.sent.at(-1))&&Object.isFrozen(c.sent.at(-1).options),'Job snapshots must be immutable');
  reports.push('Explicit method survives preset defaults; one immutable request snapshot');
  c.api.setBusy(false);
  c.api.setScan({spec:{basis:'sto-3g',count:1},options:{method:'rhf',multiplicity:1,guess:'core',charge:0,maxIterations:100,diis:true},points:[{molecule:{basis:'sto-3g',atoms:[]},result:{}}]});
  c.node('method').value='uhf';c.node('method').emit('change');c.api.run();
  const before=c.node('method').value;
  c.api.selectPoint(0);
  equal(c.node('method').value,before,'Old scan point must not rewrite an active single point');
  equal(c.sent.at(-1).options.method,'uhf','Active job remains UHF');
  reports.push('Old scan inspection cannot change a running calculation');
  c.api.setBusy(false);
  c.node('method').value='mp2';c.node('method').emit('change');
  c.node('basis').value='6-31g';c.node('basis').emit('change');
  c.node('charge').value='1';c.node('charge').emit('change');
  c.api.selectPoint(0);
  equal(c.node('method').value,'mp2','Inspecting an old point preserves a newly selected method');
  equal(c.node('basis').value,'6-31g','Inspecting an old point preserves a newly selected basis');
  equal(c.node('charge').value,'1','Inspecting an old point preserves a newly selected charge');
  c.shown.length=0;
  reports.push('Idle scan inspection preserves newer method, basis and charge edits');
  c.api.setBusy(false);
  const stale=c.sent.at(-1).id;
  c.node('method').value='mp2';c.node('method').emit('change');
  c.env.worker.onmessage({data:{type:'result',id:stale,protocol:2,buildId:'test-build',result:{method:'UHF'}}});
  equal(c.shown.length,0,'Late result after editing settings must be discarded');
  reports.push('Late result cannot repopulate edited settings');
  const reopened=controller(c.memory);
  equal(reopened.node('method').value,'mp2','Reopen restores complete selected method');
  equal(reopened.node('preset').value,'custom','Reopen restores the inspected geometry and edited setup');
  equal(reopened.node('charge').value,'1','Reopen restores charge with geometry and method');
  reopened.api.run();equal(reopened.sent.at(-1).options.method,'mp2','Reopened request uses restored MP2');
  reports.push('Reopen restores one coherent draft instead of partial browser form state');
  const analysis=controller();analysis.api.run();analysis.api.setBusy(false);
  const revisionBefore=analysis.api.snapshot().jobId;
  for(const id of ['matrix-kind','pes-energy-unit','atom-labels']) {
    analysis.node(id).emit('input');analysis.node(id).emit('change');
  }
  equal(analysis.api.snapshot().jobId,revisionBefore,'Analysis and viewer controls must not invalidate calculations');
  reports.push('Analysis and viewer controls preserve the current result');
  const menu=controller();menu.api.preview();const assignments=[];
  const native=menu.node('method');let nativeDisabled=native.disabled;
  Object.defineProperty(native,'disabled',{get:()=>nativeDisabled,set:v=>{assignments.push(v);nativeDisabled=v;}});
  native.value='uhf';native.emit('input');menu.api.flush();native.emit('change');
  equal(assignments.length,0,'Idle native select must not be disabled/re-enabled during selection');
  equal(native.value,'uhf','Committed native selection remains UHF');
  reports.push('Native input/change sequence leaves enabled select state untouched');
  const fallback=controller();fallback.node('preset').value='h2plus';fallback.node('preset').emit('change');
  equal(fallback.node('method').value,'uhf','Untouched method can adopt the ion default');
  fallback.node('method').value='mp2';fallback.node('method').emit('change');
  fallback.node('water').emit('click');equal(fallback.node('method').value,'rhf','Reset explicitly resets all defaults');
  reports.push('Initial ion defaults and explicit full Reset remain deterministic');
  const rq=settings.captureRequest(menu.form,{id:99,protocol:2,buildId:'test-build'});
  let mismatch=false;try{settings.assertResultMatches(rq,{method:'RHF',molecule:{basis:rq.basis},rhf:{charge:0}});}catch{mismatch=true;}
  check(mismatch,'Wrong-method result must be rejected');
  menu.node('method').value='rhf';equal(rq.options.method,'uhf','DOM changes cannot mutate a captured job');
  reports.push('Wrong-method result rejected; captured settings independent of form');
  return reports;
}

export async function runWorkerTests(workerSource,electronicSource) {
  const sent=[],calls=[];
  const electronic=Function('computeRHF','computeUHF','validateRHFOptions','validateUHFOptions',
    electronicSource.replace(/^import .*;$/gm,'').replace(/^export /gm,'')+'\nreturn {computeElectronic};')(
      async(_e,_m,o)=>{calls.push(o.mp2?'MP2':'RHF');return {rhf:{charge:o.charge},...(o.mp2?{mp2:{status:'complete',energy:-1.1}}:{})};},
      async(_e,_m,o)=>{calls.push('UHF');return {uhf:{charge:o.charge,multiplicity:o.multiplicity}};},
      (_m,o)=>{if(o.charge===1)throw new Error('RHF requires a positive even electron count.');return {charge:o.charge};},
      (_m,o)=>({charge:o.charge,multiplicity:o.multiplicity}));
  const self={postMessage:data=>sent.push(data)};
  let wasmURL='';
  const env={self,BUILD_ID:'test-build',WORKER_PROTOCOL:2,URL,performance:{now:()=>0},ArrayBuffer,Set,Object,
    createQChem:async config=>{wasmURL=config.locateFile('qc_core.wasm','https://test.invalid/');return {};},
    parseGeometry:(_g,_u,basis)=>({basis}),
    computeIntegrals:()=>({nao:2}),computeElectronic:electronic.computeElectronic,
    runScan:async()=>({spec:{method:'uhf'},cancelled:false,requested:1,completed:1})};
  await Function('env','with(env){return (async()=>{'+workerSource.replace(/^import .*;$/gm,'').replaceAll('import.meta.url',"'https://test.invalid/worker.mjs'")+'})()}')(env);
  equal(sent[0].protocol,2,'Worker ready handshake protocol');
  equal(sent[0].buildId,'test-build','Worker ready handshake build');
  check(wasmURL.endsWith('qc_core.wasm?v=test-build'),'WASM fetch must use the same build ID');
  let id=0;
  async function request(options,extra={}) {await self.onmessage({data:{id:++id,protocol:2,buildId:'test-build',geometry:'H 0 0 0',unit:'bohr',basis:'sto-3g',mode:'rhf',options,...extra}});return sent.at(-1);}
  for(const method of ['rhf','uhf','mp2','uhf','rhf','mp2']) {
    const response=await request({method,charge:0,multiplicity:1});
    equal(response.type,'result','Explicit method returns a result');
    equal(response.result.method,method.toUpperCase(),'Worker dispatcher must use the selected method');
    if(method==='mp2')check(!!response.result.mp2,'MP2 result must include its correlation calculation');
  }
  const count=calls.length;
  equal((await request({charge:0})).type,'error','Missing method must not fall back to RHF');
  equal((await request({method:'uhf',charge:1},{protocol:1})).type,'error','Old protocol must be rejected');
  equal((await request({method:'uhf',charge:1},{buildId:'old'})).type,'error','Mixed build must be rejected');
  equal(calls.length,count,'Invalid or stale requests must not start a solver');
  const odd=await request({method:'mp2',charge:1});equal(odd.type,'error','Odd-electron MP2 validation');check(odd.message.startsWith('MP2 requires'),'MP2 validation explains its RHF reference');
  return ['RHF/UHF/MP2 repeated worker dispatch','No silent RHF fallback for missing method','Page/worker/WASM build consistency','Explicit MP2-reference validation'];
}

if(import.meta.url===new URL(process.argv[1], 'file:').href) {
  const read=path=>fs.readFileSync(new URL(path,import.meta.url),'utf8');
  for(const report of runStateTests(read('../web/app.mjs'),read('../web/settings.mjs')))console.log(report+' PASS');
  for(const report of await runWorkerTests(read('../web/worker.mjs'),read('../web/electronic.mjs')))console.log(report+' PASS');
}
