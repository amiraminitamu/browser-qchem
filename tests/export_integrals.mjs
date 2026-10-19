// Used by the independent Python comparisons. Input is JSON on stdin.
import fs from 'node:fs';
import createQChem from '../build/wasm/dist/qc_core.mjs';
import { parseGeometry } from '../web/molecule.mjs';
import { computeIntegrals, serializable } from '../web/engine.mjs';
const input = JSON.parse(fs.readFileSync(0, 'utf8'));
const engine = await createQChem();
console.log(JSON.stringify(serializable(computeIntegrals(engine, parseGeometry(input.geometry, input.unit, input.basis)))));
