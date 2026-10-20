import fs from 'node:fs';
import createQChem from '../build/wasm/dist/qc_core.mjs';
import {parseGeometry} from '../web/molecule.mjs';
import {computeElectronic} from '../web/electronic.mjs';
import {serializable,prepareIntegrals} from '../web/engine.mjs';
import {validateUHFOptions} from '../web/uhf.mjs';
const input=JSON.parse(fs.readFileSync(0,'utf8')),engine=await createQChem();
const molecule=parseGeometry(input.geometry,input.unit,input.basis);
if(process.argv.includes('--initial-uhf')) {
  const s=validateUHFOptions(molecule,input.options),{handle,result}=prepareIntegrals(engine,molecule);
  const u=engine._qc_uhf_create(handle,s.nalpha,s.nbeta,s.maxIterations,1,s.guess==='mix'?1:0,1e-10,1e-8,1e-8);
  if(!u) throw new Error(engine.UTF8ToString(engine._qc_uhf_error()));
  try {console.log(JSON.stringify([0,1].map(k=>{const ptr=engine._qc_uhf_data(u,k);return Array.from(engine.HEAPF64.slice(ptr/8,ptr/8+result.nao**2));})));}
  finally {engine._qc_uhf_free(u);engine._qc_free(handle);}
} else console.log(JSON.stringify(serializable(await computeElectronic(engine,molecule,input.options))));
