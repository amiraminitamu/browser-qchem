// SPDX-License-Identifier: AGPL-3.0-only
import {assertResources} from './resources.mjs?v=e933e699b7cf09943ce2';
import { packBasis } from './molecule.mjs?v=e933e699b7cf09943ce2';

// Copy results before freeing C memory. Allocations may replace WASM heap views.
export function prepareIntegrals(engine, molecule, dense=true, options={}) {
  if(dense)assertResources(molecule,{...options,mode:'integrals'});
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
    handle = engine._qc_compute_mode(dense?1:0,data.charges.length, charges, xyz, data.shells.length / 3,
      shells, data.exponents.length, exponents, coefficients);
    if (!handle) throw new Error(engine.UTF8ToString(engine._qc_error()));
    const n = engine._qc_nao(handle);
    const result = { molecule, nao: n, nshell: data.shells.length / 3, aoLabels: data.aoLabels,
      nuclearRepulsion: engine._qc_enuc(handle), integralUnits: 'atomic units', layout: 'row-major', eriNotation: '(ij|kl)' };
    ['overlap', 'kinetic', 'nuclear', 'eri'].forEach((name, kind) => {
      if(kind===3&&!dense)return;
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

export function computeIntegrals(engine, molecule, options={}) {
  const resources=assertResources(molecule,{...options,mode:'integrals'});
  const { handle, result } = prepareIntegrals(engine, molecule,true,options);
  engine._qc_free(handle);
  result.resources=resources;
  return result;
}

export function checkInvariants(result) {
  const n = result.nao;
  let diagonalError = 0, matrixSymmetryError = 0, eriSymmetryError = 0;
  for (let i = 0; i < n; ++i) {
    diagonalError = Math.max(diagonalError, Math.abs(result.overlap[i*n+i] - 1));
    for (let j = 0; j < n; ++j) {
      for (const matrix of [result.overlap, result.kinetic, result.nuclear]) matrixSymmetryError = Math.max(matrixSymmetryError, Math.abs(matrix[i*n+j] - matrix[j*n+i]));
      if(result.eri) for (let k = 0; k < n; ++k) for (let l = 0; l < n; ++l) {
        const value = result.eri[((i*n+j)*n+k)*n+l];
        for (const other of [((j*n+i)*n+k)*n+l, ((i*n+j)*n+l)*n+k, ((k*n+l)*n+i)*n+j]) eriSymmetryError = Math.max(eriSymmetryError, Math.abs(value - result.eri[other]));
      }
    }
  }
  return { diagonalError, matrixSymmetryError, eriSymmetryError:result.eri?eriSymmetryError:null,
    passed: Math.max(diagonalError, matrixSymmetryError, eriSymmetryError) < 1e-10 };
}
export function serializable(value) {
  if(ArrayBuffer.isView(value)) return Array.from(value);
  if(Array.isArray(value)) return value.map(serializable);
  if(value&&typeof value==='object') return Object.fromEntries(Object.entries(value).map(([k,v])=>[k,serializable(v)]));
  return value;
}

// Correlation fitting is independent of the conventional SCF integral handle.
export function prepareDFIntegrals(engine,molecule,auxiliary) {
  const orbital=packBasis(molecule),allocations=[];
  const copy=(values,integer=false)=>{
    const pointer=engine._malloc(values.length*(integer?4:8));
    if(!pointer)throw new Error('Could not allocate DF input memory.');
    allocations.push(pointer);(integer?engine.HEAP32:engine.HEAPF64).set(values,pointer/(integer?4:8));return pointer;
  };
  try {
    const charges=copy(orbital.charges,true),xyz=copy(orbital.xyz),shells=copy(orbital.shells,true),exponents=copy(orbital.exponents),coefficients=copy(orbital.coefficients);
    const auxshells=copy(auxiliary.shells,true),auxexponents=copy(auxiliary.exponents),auxcoefficients=copy(auxiliary.coefficients);
    const handle=engine._qc_df_compute(orbital.charges.length,charges,xyz,orbital.shells.length/3,shells,orbital.exponents.length,exponents,coefficients,auxiliary.shells.length/3,auxshells,auxiliary.exponents.length,auxexponents,auxcoefficients);
    if(!handle)throw new Error(engine.UTF8ToString(engine._qc_df_error()));
    return handle;
  }finally{allocations.forEach(pointer=>engine._free(pointer));}
}

// Exact ERI inspector is independent of the SCF backend and allocates only S/T/V.
export function computeERI(engine,molecule,indices) {
  if(!Array.isArray(indices)||indices.length!==4)throw new Error('Choose four AO indices.');
  const {handle,result}=prepareIntegrals(engine,molecule,false);
  try {
    if(indices.some(i=>!Number.isInteger(i)||i<0||i>=result.nao))throw new Error('AO index is outside this basis.');
    const value=engine._qc_eri(handle,...indices);
    if(!Number.isFinite(value))throw new Error('ERI evaluation failed.');
    return value;
  }finally{engine._qc_free(handle);}
}
