// SPDX-License-Identifier: AGPL-3.0-only
import {METHOD_SPECS} from './methods.mjs?v=e933e699b7cf09943ce2';
import createQChem from './qc_core.mjs?v=e933e699b7cf09943ce2';
import { parseGeometry } from './molecule.mjs?v=e933e699b7cf09943ce2';
import { computeIntegrals,computeERI } from './engine.mjs?v=e933e699b7cf09943ce2';
import { computeElectronic } from './electronic.mjs?v=e933e699b7cf09943ce2';
import { runScan } from './scan.mjs?v=e933e699b7cf09943ce2';
import {BUILD_ID,WORKER_PROTOCOL} from './build-info.mjs?v=e933e699b7cf09943ce2';
let active=null;
const transfers=result=>[...new Set([...Object.values(result),...Object.values(result.rhf??{}),...Object.values(result.uhf??{})]
  .filter(ArrayBuffer.isView).map(value=>value.buffer))];
try {
  const engine=await createQChem({locateFile:(path,prefix)=>{
    const url=new URL(path,prefix||import.meta.url);url.searchParams.set('v',BUILD_ID);return url.href;
  }});
  self.onmessage=async ({data})=>{
    if (data.type==='cancel') { if (active&&active.id===data.id) active.cancelled=true; return; }
    if(data.type==='eri-query') {
      const reply={type:'eri-value',query:data.query,protocol:WORKER_PROTOCOL,buildId:BUILD_ID};
      try {
        if(data.protocol!==WORKER_PROTOCOL||data.buildId!==BUILD_ID)throw new Error('Page and engine versions differ.');
        if(active)throw new Error('Wait for the calculation to finish before inspecting ERIs.');
        const molecule=parseGeometry(data.geometry,data.unit,data.basis);
        self.postMessage({...reply,value:computeERI(engine,molecule,data.indices)});
      }catch(error){self.postMessage({...reply,error:error.message});}
      return;
    }
    if (active) { self.postMessage({type:'error',id:data.id,protocol:WORKER_PROTOCOL,buildId:BUILD_ID,message:'A calculation is already running.'}); return; }
    const job={id:data.id??0,cancelled:false}; active=job;
    const send=(type,payload={},transfer=[])=>self.postMessage({type,id:job.id,protocol:WORKER_PROTOCOL,buildId:BUILD_ID,...payload},transfer);
    try {
      if(data.protocol!==WORKER_PROTOCOL||data.buildId!==BUILD_ID)throw new Error('Page and engine versions differ. Reload the page.');
      if(!['rhf','scan','integrals'].includes(data.mode))throw new Error('Calculation type is missing or invalid.');
      if(data.mode!=='integrals'&&!Object.hasOwn(METHOD_SPECS,data.options?.method))throw new Error('Electronic method is missing or invalid. Choose a supported electronic method and reload if necessary.');
      const start=performance.now();
      const molecule=parseGeometry(data.geometry,data.unit,data.basis);
      const options={...data.options,isCancelled:()=>job.cancelled};
      if (data.mode==='scan') {
        const scan=await runScan(engine,molecule,data.scan,options,
          (point,spec)=>send('point',{point,spec}),
          ({index,entry})=>send('iteration',{index,entry}),
          progress=>send('post-scf',progress));
        send('scan-complete',{spec:scan.spec,summary:{cancelled:scan.cancelled,requested:scan.requested,completed:scan.completed,elapsedMs:performance.now()-start}});
      } else {
        const result=data.mode==='integrals'?computeIntegrals(engine,molecule,options)
          :await computeElectronic(engine,molecule,options,entry=>send('iteration',{entry}),progress=>send('post-scf',{progress}));
        result.elapsedMs=performance.now()-start;
        if (job.cancelled) send('cancelled');
        else send('result',{result},transfers(result));
      }
    } catch(error) {
      send(error.name==='AbortError'?'cancelled':'error',{message:error.message});
    } finally { active=null; }
  };
  self.postMessage({type:'ready',protocol:WORKER_PROTOCOL,buildId:BUILD_ID});
} catch(error) { self.postMessage({type:'fatal',protocol:WORKER_PROTOCOL,buildId:BUILD_ID,message:error.message}); }

