import fs from 'node:fs';
import createQChem from '../build/wasm/dist/qc_core.mjs';
import {parseGeometry} from '../web/molecule.mjs';
import {runScan} from '../web/scan.mjs';
import {serializable} from '../web/engine.mjs';
const input=JSON.parse(fs.readFileSync(0,'utf8')),engine=await createQChem();
input.options={scfIntegrals:'exact',...input.options};
console.log(JSON.stringify(serializable(await runScan(engine,parseGeometry(input.geometry,input.unit,input.basis),input.scan,input.options))));
