"""Regenerate the bundled H-Ne catalog with PySCF 2.14.0 (development only)."""
import json
from pathlib import Path
import pyscf
from pyscf import gto

if pyscf.__version__ != '2.14.0':
    raise RuntimeError('Use PySCF 2.14.0 to reproduce the pinned basis tables.')
catalog = {}
for name,label in [('sto-3g','STO-3G'),('6-31g','6-31G'),('cc-pvdz','cc-pVDZ')]:
    elements = {}
    for charge,symbol in enumerate(['H','He','Li','Be','B','C','N','O','F','Ne'],1):
        shells = []
        for shell in gto.basis.load(name,symbol):
            angular,primitives = shell[0],shell[1:]
            for contraction in range(1,len(primitives[0])):
                shells.append({'l':angular,'exponents':[p[0] for p in primitives],
                               'coefficients':[p[contraction] for p in primitives]})
        elements[symbol] = {'charge':charge,'shells':shells}
    catalog[name] = {'label':label,'elements':elements}
source = '// Scientific basis parameters transcribed from PySCF 2.14.0.\n'
source += '// General contractions are split into single-contraction shells in the same AO order.\n'
source += '// See README for scientific references and data provenance.\n'
source += 'export const BASIS_CATALOG = '+json.dumps(catalog,separators=(',',':'))+';\n'
source += 'export const STO3G = BASIS_CATALOG["sto-3g"].elements;\n'
source += "export {MAX_AO} from './limits.mjs';\n"
source += 'export function basisLabel(name) { return BASIS_CATALOG[name]?.label ?? name; }\n'
(Path(__file__).resolve().parents[1]/'web/basis.mjs').write_text(source)
