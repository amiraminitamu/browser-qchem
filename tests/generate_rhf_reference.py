"""Generate a water RHF reference using PySCF only."""
import json
from pathlib import Path
import pyscf
from pyscf import gto,scf
mol=gto.M(atom='O 0 0 0; H 0 -.757 .587; H 0 .757 .587',basis='sto-3g',unit='Angstrom',verbose=0)
mf=scf.RHF(mol)
mf.conv_tol=1e-12; mf.conv_tol_grad=1e-10; mf.init_guess='1e'; mf.kernel()
assert mf.converged
reference={'source':f'PySCF {pyscf.__version__}','energy':mf.e_tot,'nao':mol.nao_nr(),
    'density':mf.make_rdm1().reshape(-1).tolist(),'orbitalEnergies':mf.mo_energy.tolist(),
    'fock':mf.get_fock().reshape(-1).tolist()}
path=Path(__file__).resolve().parent/'fixtures/water-rhf-sto3g.json'
path.write_text(json.dumps(reference,indent=2)+'\n')
print('Water reference:',mf.e_tot)
