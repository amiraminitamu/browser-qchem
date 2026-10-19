"""Compare every S, T, V and ERI entry with PySCF, natively or through WASM.

python3 tests/check_molecular.py --native build/native/libqc_core.dylib
python3 tests/check_molecular.py --wasm
"""
import argparse
import ctypes as ct
import json
from pathlib import Path
import subprocess

import numpy as np
from pyscf import gto

ROOT = Path(__file__).resolve().parents[1]
NAMES = ('overlap', 'kinetic', 'nuclear', 'eri')
OPERATORS = ('int1e_ovlp_sph', 'int1e_kin_sph', 'int1e_nuc_sph', 'int2e_sph')

def native_result(library, mol):
    charges = np.ascontiguousarray(mol.atom_charges(), dtype=np.int32)
    xyz = np.ascontiguousarray(mol.atom_coords(), dtype=np.float64)
    shells, exponents, coefficients = [], [], []
    for atom in range(mol.natm):
        for shell in mol._basis[mol.atom_symbol(atom)]:
            l, primitives = shell[0], shell[1:]
            for contraction in range(1, len(primitives[0])):
                shells.extend([atom, l, len(primitives)])
                exponents.extend(p[0] for p in primitives)
                coefficients.extend(p[contraction] for p in primitives)
    shells = np.asarray(shells, dtype=np.int32)
    exponents = np.asarray(exponents, dtype=np.float64)
    coefficients = np.asarray(coefficients, dtype=np.float64)
    pointer = library.qc_compute(mol.natm, charges.ctypes.data_as(ct.POINTER(ct.c_int)),
        xyz.ctypes.data_as(ct.POINTER(ct.c_double)), len(shells)//3,
        shells.ctypes.data_as(ct.POINTER(ct.c_int)), len(exponents),
        exponents.ctypes.data_as(ct.POINTER(ct.c_double)), coefficients.ctypes.data_as(ct.POINTER(ct.c_double)))
    if not pointer:
        raise RuntimeError(library.qc_error().decode())
    try:
        n = library.qc_nao(pointer)
        result = {'nao': n, 'nuclearRepulsion': library.qc_enuc(pointer)}
        for kind, name in enumerate(NAMES):
            result[name] = np.ctypeslib.as_array(library.qc_data(pointer, kind), shape=(n**(4 if kind == 3 else 2),)).copy()
        return result
    finally:
        library.qc_free(pointer)

def load_library(path):
    library = ct.CDLL(str(Path(path).resolve()))
    ip, dp = ct.POINTER(ct.c_int), ct.POINTER(ct.c_double)
    library.qc_compute.argtypes = [ct.c_int,ip,dp,ct.c_int,ip,ct.c_int,dp,dp]
    library.qc_compute.restype = ct.c_void_p
    library.qc_error.restype = ct.c_char_p
    library.qc_nao.argtypes = [ct.c_void_p]; library.qc_nao.restype = ct.c_int
    library.qc_enuc.argtypes = [ct.c_void_p]; library.qc_enuc.restype = ct.c_double
    library.qc_data.argtypes = [ct.c_void_p,ct.c_int]; library.qc_data.restype = dp
    library.qc_free.argtypes = [ct.c_void_p]; library.qc_free.restype = None
    return library

def cases():
    for symbol in ['H','He','Li','Be','B','C','N','O','F','Ne']:
        yield symbol, f'{symbol} 0 0 0', 'Bohr'
    yield 'H2', 'H 0 0 0; H 0 0 1.4', 'Bohr'
    yield 'water', 'O 0 0 0; H 0 -0.757 0.587; H 0 0.757 0.587', 'Angstrom'
    # A skew geometry exercises all p directions and nonzero off-diagonal blocks.
    yield 'skew water', 'O 0.14 -0.23 0.31; H 0.79 0.02 0.91; H -0.51 0.61 0.44', 'Angstrom'
    yield 'methane', 'C 0 0 0; H .63 .63 .63; H -.63 -.63 .63; H -.63 .63 -.63; H .63 -.63 -.63', 'Angstrom'
    yield 'ammonia', 'N 0 0 .1; H .93 0 -.25; H -.465 .8054 -.25; H -.465 -.8054 -.25', 'Angstrom'
    yield 'HF', 'H 0 0 0; F .2 .1 .9', 'Angstrom'
    yield 'LiH', 'Li 0 0 0; H 0 0 1.6', 'Angstrom'
    yield 'BeH2', 'Be 0 0 0; H 0 0 1.3; H 0 0 -1.3', 'Angstrom'

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--native'); mode.add_argument('--wasm', action='store_true')
    args = parser.parse_args()
    library = load_library(args.native) if args.native else None
    largest = 0.0
    for name, geometry, unit in cases():
        # Odd-electron atoms need a spin value for PySCF molecule construction.
        electron_count = sum(gto.charge(line.split()[0]) for line in geometry.split(';'))
        mol = gto.M(atom=geometry, unit=unit, basis='sto-3g', spin=electron_count % 2, verbose=0)
        if library:
            actual = native_result(library, mol)
        else:
            process = subprocess.run(['node', str(ROOT/'tests/export_integrals.mjs')],
                input=json.dumps({'geometry':geometry.replace(';','\n'), 'unit':unit.lower()}),
                text=True, capture_output=True, check=True, cwd=ROOT)
            actual = json.loads(process.stdout)
        assert actual['nao'] == mol.nao_nr()
        errors = []
        for key, operator in zip(NAMES, OPERATORS):
            reference = mol.intor(operator, aosym='s1').reshape(-1)
            values = np.asarray(actual[key]).reshape(-1)
            assert values.shape == reference.shape, (name,key,values.shape,reference.shape)
            assert np.all(np.isfinite(values))
            error = float(np.max(np.abs(values-reference)))
            assert error < 1e-10, (name,key,error)
            errors.append(error)
        enuc_error = abs(actual['nuclearRepulsion'] - mol.energy_nuc())
        assert enuc_error < 1e-11, (name,enuc_error)
        largest = max(largest,*errors)
        print(f'{name:12s} {mol.nao_nr():2d} AOs | S {errors[0]:.2e} T {errors[1]:.2e} V {errors[2]:.2e} ERI {errors[3]:.2e} PASS')
    if library:
        # A custom d shell tests the generic C interface beyond the UI's STO-3G.
        mol = gto.M(atom='H 0.2 -0.3 0.4; H -0.4 0.7 1.1',unit='Bohr',
            basis={'H':[[0,[1.1,1.0]],[1,[.7,.8],[.2,.3]],[2,[.6,1.0]]]},verbose=0)
        actual = native_result(library,mol)
        for name,op in zip(NAMES,OPERATORS):
            error = float(np.max(np.abs(actual[name] - mol.intor(op,aosym='s1').reshape(-1))))
            assert error < 1e-10, (name,error)
            largest=max(largest,error)
        print('Custom s/p/d shells: PASS')
    print(f'Maximum absolute integral error: {largest:.3e}')

if __name__ == '__main__':
    main()
