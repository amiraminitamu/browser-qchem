// SPDX-License-Identifier: AGPL-3.0-only
import { packBasis } from './molecule.mjs';

// Copy results before freeing C memory. Allocations may replace WASM heap views.
export function prepareIntegrals(engine, molecule) {
  const data = packBasis(molecule), allocations = [];
  let handle = 0;
  function copyIn(values, integer = false) {
    const pointer = engine._malloc(values.length * (integer ? 4 : 8));
    if (!pointer) throw new Error('Could not allocate WebAssembly input memory.');
    allocations.push(pointer);
    if (integer) engine.HEAP32.set(values, pointer / 4);
    else engine.HEAPF64.set(values, pointer / 8);
    return pointer;
  }
  try {
    const charges = copyIn(data.charges, true), xyz = copyIn(data.xyz);
    const shells = copyIn(data.shells, true), exponents = copyIn(data.exponents), coefficients = copyIn(data.coefficients);
    handle = engine._qc_compute(data.charges.length, charges, xyz, data.shells.length / 3,
      shells, data.exponents.length, exponents, coefficients);
    if (!handle) throw new Error(engine.UTF8ToString(engine._qc_error()));
    const n = engine._qc_nao(handle);
    const result = { molecule, nao: n, nshell: data.shells.length / 3, aoLabels: data.aoLabels,
      nuclearRepulsion: engine._qc_enuc(handle), integralUnits: 'atomic units', layout: 'row-major', eriNotation: '(ij|kl)' };
    ['overlap', 'kinetic', 'nuclear', 'eri'].forEach((name, kind) => {
      const pointer = engine._qc_data(handle, kind), size = kind === 3 ? n ** 4 : n ** 2;
      result[name] = engine.HEAPF64.slice(pointer / 8, pointer / 8 + size);
      if (!result[name].every(Number.isFinite)) throw new Error(`Nonfinite ${name} integral.`);
    });
    result.checks = checkInvariants(result);
    return { handle, result };
  } catch (error) {
    if (handle) engine._qc_free(handle);
    throw error;
  } finally {
    allocations.forEach(pointer => engine._free(pointer));
  }
}

export function computeIntegrals(engine, molecule) {
  const { handle, result } = prepareIntegrals(engine, molecule);
  engine._qc_free(handle);
  return result;
}

export function checkInvariants(result) {
  const n = result.nao;
  let diagonalError = 0, matrixSymmetryError = 0, eriSymmetryError = 0;
  for (let i = 0; i < n; ++i) {
    diagonalError = Math.max(diagonalError, Math.abs(result.overlap[i*n+i] - 1));
    for (let j = 0; j < n; ++j) {
      for (const matrix of [result.overlap, result.kinetic, result.nuclear]) matrixSymmetryError = Math.max(matrixSymmetryError, Math.abs(matrix[i*n+j] - matrix[j*n+i]));
      for (let k = 0; k < n; ++k) for (let l = 0; l < n; ++l) {
        const value = result.eri[((i*n+j)*n+k)*n+l];
        for (const other of [((j*n+i)*n+k)*n+l, ((i*n+j)*n+l)*n+k, ((k*n+l)*n+i)*n+j]) eriSymmetryError = Math.max(eriSymmetryError, Math.abs(value - result.eri[other]));
      }
    }
  }
  return { diagonalError, matrixSymmetryError, eriSymmetryError,
    passed: Math.max(diagonalError, matrixSymmetryError, eriSymmetryError) < 1e-10 };
}
export function serializable(value) {
  if(ArrayBuffer.isView(value)) return Array.from(value);
  if(Array.isArray(value)) return value.map(serializable);
  if(value&&typeof value==='object') return Object.fromEntries(Object.entries(value).map(([k,v])=>[k,serializable(v)]));
  return value;
}
