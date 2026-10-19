// SPDX-License-Identifier: AGPL-3.0-only
import createQChem from './qc_core.mjs';
import { parseGeometry } from './molecule.mjs';
import { computeIntegrals } from './engine.mjs';
import { computeRHF } from './rhf.mjs';
import { runScan } from './scan.mjs';
let active=null;
const transfers=result=>[...new Set([...Object.values(result),...Object.values(result.rhf??{})]
  .filter(ArrayBuffer.isView).map(value=>value.buffer))];
try {
  const engine=await createQChem();
  self.onmessage=async ({data})=>{
    if (data.type==='cancel') { if (active&&active.id===data.id) active.cancelled=true; return; }
    if (active) { self.postMessage({type:'error',id:data.id,message:'A calculation is already running.'}); return; }
    const job={id:data.id??0,cancelled:false}; active=job;
    const send=(type,payload={})=>self.postMessage({type,id:job.id,...payload});
    try {
      const start=performance.now();
      const molecule=parseGeometry(data.geometry,data.unit,data.basis);
      const options={...data.options,isCancelled:()=>job.cancelled};
      if (data.mode==='scan') {
        const scan=await runScan(engine,molecule,data.scan,options,
          (point,spec)=>send('point',{point,spec}),
          ({index,entry})=>send('iteration',{index,entry}));
        send('scan-complete',{spec:scan.spec,summary:{cancelled:scan.cancelled,requested:scan.requested,completed:scan.completed,elapsedMs:performance.now()-start}});
      } else {
        const result=data.mode==='integrals'?computeIntegrals(engine,molecule)
          :await computeRHF(engine,molecule,options,entry=>send('iteration',{entry}));
        result.elapsedMs=performance.now()-start;
        if (job.cancelled) send('cancelled');
        else self.postMessage({type:'result',id:job.id,result},transfers(result));
      }
    } catch(error) {
      send(error.name==='AbortError'?'cancelled':'error',{message:error.message});
    } finally { active=null; }
  };
  self.postMessage({type:'ready'});
} catch(error) { self.postMessage({type:'fatal',message:error.message}); }
