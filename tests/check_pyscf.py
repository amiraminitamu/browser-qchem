"""Optional independent environment check: native wrapper versus PySCF.

Usage: python3 tests/check_pyscf.py build/native/libqc_core.dylib  # macOS
       python3 tests/check_pyscf.py build/native/libqc_core.so     # Linux
Requires PySCF only for this development check, never for the webpage.
"""
import ctypes
import math
import sys
from pathlib import Path
from pyscf import gto

library = ctypes.CDLL(str(Path(sys.argv[1]).resolve()))
library.h2_overlap.argtypes = [ctypes.c_double]
library.h2_overlap.restype = ctypes.c_double
for r in [0, 0.1, 0.7, 1.4, 2, 4, 10, 20]:
    mol = gto.M(atom=[('H', (0, 0, 0)), ('H', (0, 0, r))],
                basis='sto-3g', unit='Bohr', verbose=0)
    reference = mol.intor('int1e_ovlp_sph')[0, 1]
    actual = library.h2_overlap(r)
    if not math.isfinite(actual) or abs(actual - reference) >= 1e-12:
        raise AssertionError(f'{r} bohr: {actual} != {reference}')
    print(f'{r:4.1f} bohr: {actual:.15f}, PySCF error {abs(actual-reference):.3e} PASS')
