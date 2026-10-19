import fs from 'node:fs';
import createQChem from '../build/wasm/dist/qc_core.mjs';
import { parseGeometry } from '../web/molecule.mjs';
import { computeRHF } from '../web/rhf.mjs';
import { serializable } from '../web/engine.mjs';
const input = JSON.parse(fs.readFileSync(0,'utf8'));
const engine = await createQChem();
console.log(JSON.stringify(serializable(await computeRHF(engine,parseGeometry(input.geometry,input.unit,input.basis),input.options))));
