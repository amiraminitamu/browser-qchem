"""Bundle def2-universal-JKFIT H-Ne, independent of the correlation RI basis."""
import json
from pathlib import Path
import pyscf
from pyscf import gto
assert pyscf.__version__ == '2.14.0'
elements={}
for symbol in ['H','He','Li','Be','B','C','N','O','F','Ne']:
    shells=[]
    for angular,*primitives in gto.basis.load('def2-universal-jkfit',symbol):
        for c in range(1,len(primitives[0])):
            shells.append({'l':angular,'exponents':[p[0] for p in primitives],'coefficients':[p[c] for p in primitives]})
    elements[symbol]={'shells':shells}
p=Path(__file__).resolve().parents[1]/'web/jkfit-basis.mjs'
p.write_text('// SPDX-License-Identifier: AGPL-3.0-only\n// def2-universal-JKFIT parameters from PySCF 2.14.0; spherical H-Ne.\nexport const JKFIT_BASIS='+json.dumps({'id':'def2-universal-jkfit','label':'def2-universal-JKFIT','elements':elements},separators=(',',':'))+';\n')
