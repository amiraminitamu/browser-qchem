"""Validate the native UHF/MP2 C APIs against the committed PySCF fixtures.
python3 tests/check_methods_native.py build/native/libqc_core.so
"""
import ctypes as ct,json,sys
from pathlib import Path
from pyscf import gto,lib
from check_molecular import load_library,native_handle
from check_rhf import configure
lib.num_threads(1)
library=configure(load_library(sys.argv[1]));ROOT=Path(__file__).resolve().parents[1]
library.qc_uhf_create.argtypes=[ct.c_void_p,*([ct.c_int]*5),*([ct.c_double]*3)];library.qc_uhf_create.restype=ct.c_void_p
library.qc_mp2_create.argtypes=[ct.c_void_p];library.qc_mp2_create.restype=ct.c_void_p
for prefix in ['uhf','mp2']:
    for suffix in ['step','status']:
        f=getattr(library,f'qc_{prefix}_{suffix}');f.argtypes=[ct.c_void_p];f.restype=ct.c_int
    f=getattr(library,f'qc_{prefix}_scalar');f.argtypes=[ct.c_void_p,ct.c_int];f.restype=ct.c_double
    f=getattr(library,f'qc_{prefix}_free');f.argtypes=[ct.c_void_p];f.restype=None
    getattr(library,f'qc_{prefix}_error').restype=ct.c_char_p
cases=json.loads((ROOT/'tests/fixtures/methods.json').read_text())['cases']
for case in cases:
    options=case['options'];mol=gto.M(atom=case['geometry'],unit=case['unit'],basis=case['basis'],charge=options['charge'],spin=options['multiplicity']-1,verbose=0)
    integrals=native_handle(library,mol);solver=None;post=None;method=options['method']
    try:
        if method=='uhf':solver=library.qc_uhf_create(integrals,*mol.nelec,500,1,int(options['guess']=='mix'),1e-10,1e-8,1e-8)
        else:solver=library.qc_rhf_create(integrals,mol.nelectron,500,1,1e-10,1e-8,1e-8)
        prefix='uhf' if method=='uhf' else 'rhf';assert solver,getattr(library,f'qc_{prefix}_error')().decode()
        step=getattr(library,f'qc_{prefix}_step');status=0
        while status==0:status=step(solver)
        assert status==1,(case['name'],method,status)
        actual=getattr(library,f'qc_{prefix}_scalar')(solver,0);assert abs(actual-case['referenceEnergy'])<2e-8
        if method=='uhf':assert abs(library.qc_uhf_scalar(solver,5)-case['spinSquare'])<2e-6
        else:
            post=library.qc_mp2_create(solver);assert post;status=library.qc_mp2_status(post)
            while status==0:status=library.qc_mp2_step(post)
            assert status==1,library.qc_mp2_error().decode()
            for k,key in enumerate(['energy','correlationEnergy','oppositeSpinEnergy','sameSpinEnergy']):assert abs(library.qc_mp2_scalar(post,k)-case[key])<2e-8,(case['name'],key)
    finally:
        if post:library.qc_mp2_free(post)
        if solver:getattr(library,f'qc_{prefix}_free')(solver)
        library.qc_free(integrals)
print(f'{len(cases)} native UHF/MP2 C API comparisons PASS')
