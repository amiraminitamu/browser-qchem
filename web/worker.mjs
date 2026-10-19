// SPDX-License-Identifier: AGPL-3.0-only
import createQChem from './qc_core.mjs';
import { parseGeometry } from './molecule.mjs';
import { computeIntegrals } from './engine.mjs';
try {
  const engine = await createQChem();
  self.onmessage = ({ data }) => {
    try {
      const start = performance.now();
      const result = computeIntegrals(engine, parseGeometry(data.geometry, data.unit));
      result.elapsedMs = performance.now() - start;
      self.postMessage({ type: 'result', result },
        [result.overlap.buffer, result.kinetic.buffer, result.nuclear.buffer, result.eri.buffer]);
    } catch (error) { self.postMessage({ type: 'error', message: error.message }); }
  };
  self.postMessage({ type: 'ready' });
} catch (error) { self.postMessage({ type: 'fatal', message: error.message }); }
