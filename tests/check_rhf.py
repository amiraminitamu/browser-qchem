"""Validate RHF energy, density, Fock, orbital energies, and SCF identities.

python3 tests/check_rhf.py --native build/native/libqc_core.dylib
python3 tests/check_rhf.py --wasm
"""
import argparse
import ctypes as ct
import json
from pathlib import Path
import subprocess
import numpy as np
from pyscf import gto,scf
from check_molecular import load_library,native_handle

ROOT=Path(__file__).resolve().parents[1]

def native_rhf(library,mol,diis=True,maxiter=100):
    integrals=native_handle(library,mol)
    pointer=None
    try:
        pointer=library.qc_rhf_create(integrals,mol.nelectron,maxiter,int(diis),1e-10,1e-8,1e-8)
        if not pointer:
            raise RuntimeError(library.qc_rhf_error().decode())
        status=0
        while status==0:
            status=library.qc_rhf_step(pointer)
        if status<0:
            raise RuntimeError(library.qc_rhf_error().decode())
        n=mol.nao_nr()
        r={'converged':status==1,'energy':library.qc_rhf_scalar(pointer,0),'iterations':library.qc_rhf_iteration(pointer)}
        for kind,key in enumerate(['density','fock','coefficients','orbitalEnergies','occupations']):
            r[key]=np.ctypeslib.as_array(library.qc_rhf_data(pointer,kind),shape=(n*n if kind<3 else n,)).copy()
        return r
    finally:
        if pointer: library.qc_rhf_free(pointer)
        library.qc_free(integrals)

def configure(library):
    library.qc_rhf_create.argtypes=[ct.c_void_p,ct.c_int,ct.c_int,ct.c_int,ct.c_double,ct.c_double,ct.c_double]
    library.qc_rhf_create.restype=ct.c_void_p
    library.qc_rhf_error.restype=ct.c_char_p
    for name in ['qc_rhf_step','qc_rhf_iteration']:
        getattr(library,name).argtypes=[ct.c_void_p]; getattr(library,name).restype=ct.c_int
    library.qc_rhf_scalar.argtypes=[ct.c_void_p,ct.c_int]; library.qc_rhf_scalar.restype=ct.c_double
    library.qc_rhf_data.argtypes=[ct.c_void_p,ct.c_int]; library.qc_rhf_data.restype=ct.POINTER(ct.c_double)
    library.qc_rhf_free.argtypes=[ct.c_void_p]; library.qc_rhf_free.restype=None
    return library

CASES=[
 ('H2','H 0 0 0; H 0 0 1.4','Bohr',0),
 ('stretched H2','H 0 0 0; H 0 0 3','Bohr',0),
 ('water','O 0 0 0; H 0 -.757 .587; H 0 .757 .587','Angstrom',0),
 ('skew water','O .14 -.23 .31; H .79 .02 .91; H -.51 .61 .44','Angstrom',0),
 ('methane','C 0 0 0; H .63 .63 .63; H -.63 -.63 .63; H -.63 .63 -.63; H .63 -.63 -.63','Angstrom',0),
 ('ammonia','N 0 0 .1; H .93 0 -.25; H -.465 .8054 -.25; H -.465 -.8054 -.25','Angstrom',0),
 ('HF','H 0 0 0; F .2 .1 .9','Angstrom',0),
 ('LiH','Li 0 0 0; H 0 0 1.6','Angstrom',0),
 ('BeH2','Be 0 0 0; H 0 0 1.3; H 0 0 -1.3','Angstrom',0),
 ('He','He 0 0 0','Bohr',0),
 ('Ne','Ne 0 0 0','Bohr',0),
 ('HeH+','He 0 0 0; H 0 0 1.5','Bohr',1),
 ('H3+','H 0 0 0; H 1.65 0 0; H .825 1.428941916 0','Bohr',1),
]

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    mode=parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--native'); mode.add_argument('--wasm',action='store_true')
    args=parser.parse_args()
    library=configure(load_library(args.native)) if args.native else None
    max_energy_error=max_density_error=0.0
    for name,geometry,unit,charge in CASES:
        mol=gto.M(atom=geometry,basis='sto-3g',unit=unit,charge=charge,verbose=0)
        mol.incore_anyway=True # Small references; no OS memory probe is needed.
        ref=scf.RHF(mol)
        ref.conv_tol=1e-12; ref.conv_tol_grad=1e-10; ref.max_cycle=200; ref.init_guess='1e'
        ref.kernel()
        assert ref.converged,(name,'PySCF did not converge')
        if library:
            actual=native_rhf(library,mol)
        else:
            process=subprocess.run(['node',str(ROOT/'tests/export_rhf.mjs')],
                input=json.dumps({'geometry':geometry.replace(';','\n'),'unit':unit.lower(),'options':{'charge':charge}}),
                text=True,capture_output=True,check=True,cwd=ROOT)
            actual=json.loads(process.stdout)['rhf']
        assert actual['converged'],(name,'RHF did not converge')
        n=mol.nao_nr(); s=mol.intor('int1e_ovlp')
        d=np.asarray(actual['density']).reshape(n,n); f=np.asarray(actual['fock']).reshape(n,n)
        c=np.asarray(actual['coefficients']).reshape(n,n); eps=np.asarray(actual['orbitalEnergies'])
        occ=np.asarray(actual['occupations'])
        energy_error=abs(actual['energy']-ref.e_tot)
        density_error=float(np.max(np.abs(d-ref.make_rdm1())))
        orbital_error=float(np.max(np.abs(eps-ref.mo_energy)))
        assert energy_error<1e-9,(name,'energy',energy_error)
        assert density_error<1e-7,(name,'density',density_error)
        assert orbital_error<1e-7,(name,'orbital energies',orbital_error)
        assert np.max(np.abs(f-ref.get_fock(dm=d)))<1e-10,(name,'physical Fock')
        assert np.max(np.abs(c.T@s@c-np.eye(n)))<1e-10,(name,'orthonormality')
        assert abs(np.einsum('ij,ji',d,s)-mol.nelectron)<1e-9,(name,'electron count')
        assert np.max(np.abs(d@s@d-2*d))<1e-9,(name,'idempotency')
        assert np.max(np.abs((c*occ)@c.T-d))<1e-9,(name,'orbital density')
        assert np.max(np.abs(f@c-s@c@np.diag(eps)))<1e-7,(name,'canonical orbitals')
        assert abs(ref.energy_tot(dm=d)-actual['energy'])<1e-10,(name,'energy recomputation')
        max_energy_error=max(max_energy_error,energy_error); max_density_error=max(max_density_error,density_error)
        print(f'{name:14s} E={actual["energy"]:.12f} cycles={actual["iterations"]:2d} | dE {energy_error:.2e} dD {density_error:.2e} deps {orbital_error:.2e} PASS')
    if library:
        mol=gto.M(atom=CASES[2][1],unit='Angstrom',basis='sto-3g',verbose=0)
        assert native_rhf(library,mol,diis=False)['converged']
        assert not native_rhf(library,mol,maxiter=1)['converged']
        print('DIIS off and nonconvergence status: PASS')
    print(f'Max energy error {max_energy_error:.3e} Eh; max density entry error {max_density_error:.3e}')

if __name__=='__main__': main()
