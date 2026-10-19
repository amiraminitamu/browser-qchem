// SPDX-License-Identifier: AGPL-3.0-only
import createQChem from './qc_core.mjs';
import { parseGeometry } from './molecule.mjs';
import { computeIntegrals } from './engine.mjs';
import { computeRHF } from './rhf.mjs';
try {
  const engine = await createQChem();
  self.onmessage = async ({ data }) => {
    try {
      const start = performance.now();
      const molecule = parseGeometry(data.geometry,data.unit);
      const result = data.mode === 'rhf'
        ? await computeRHF(engine,molecule,data.options,entry => self.postMessage({ type:'iteration',entry }))
        : computeIntegrals(engine,molecule);
      result.elapsedMs = performance.now() - start;
      const buffers = [result.overlap.buffer, result.kinetic.buffer, result.nuclear.buffer, result.eri.buffer];
      if (result.rhf) for (const key of ['density','fock','coefficients','orbitalEnergies','occupations']) buffers.push(result.rhf[key].buffer);
      self.postMessage({ type:'result',result },buffers);
    } catch (error) { self.postMessage({ type: 'error', message: error.message }); }
  };
  self.postMessage({ type: 'ready' });
} catch (error) { self.postMessage({ type: 'fatal', message: error.message }); }
