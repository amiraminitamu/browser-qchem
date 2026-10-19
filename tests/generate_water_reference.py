"""Regenerate the complete water/STO-3G reference with PySCF, not our wrapper."""
import json
from pathlib import Path
import pyscf
from pyscf import gto

geometry = 'O 0 0 0; H 0 -0.757 0.587; H 0 0.757 0.587'
mol = gto.M(atom=geometry, unit='Angstrom', basis='sto-3g', verbose=0)
reference = {'source': f'PySCF {pyscf.__version__}', 'basis':'sto-3g', 'geometry':geometry,
    'unit':'angstrom', 'nao':mol.nao_nr(), 'aoLabels':mol.ao_labels(),
    'nuclearRepulsion':mol.energy_nuc(), 'layout':'row-major', 'eriNotation':'(ij|kl)'}
for name,op in [('overlap','int1e_ovlp_sph'),('kinetic','int1e_kin_sph'),('nuclear','int1e_nuc_sph'),('eri','int2e_sph')]:
    reference[name] = mol.intor(op,aosym='s1').reshape(-1).tolist()
path = Path(__file__).resolve().parent/'fixtures/water-sto3g.json'
path.parent.mkdir(exist_ok=True)
path.write_text(json.dumps(reference,indent=2)+'\n')
print(f'Wrote {path} ({mol.nao_nr()**4} ERIs)')
