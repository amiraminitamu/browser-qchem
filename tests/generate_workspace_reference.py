"""Generate small basis/scan regression references using PySCF, without our wrapper."""
import json
from pathlib import Path
import pyscf
from pyscf import gto,scf,lib
lib.num_threads(1)
cases=[]
for basis in ['sto-3g','6-31g','cc-pvdz']:
    for name,geometry in [('water','O 0 0 0;H 0 -.757 .587;H 0 .757 .587')]+[
        ('h2-'+str(r),f'H 0 0 0;H 0 0 {r}') for r in [.5,.75,1.,1.25,1.5]]:
        mol=gto.M(atom=geometry,basis=basis,unit='Angstrom',verbose=0)
        mol.incore_anyway=True
        mf=scf.RHF(mol);mf.init_guess='1e';mf.conv_tol=1e-12;mf.conv_tol_grad=1e-10;mf.max_cycle=200
        mf.kernel();assert mf.converged
        cases.append({'name':name,'basis':basis,'geometry':geometry.replace(';','\n'),
                      'unit':'angstrom','nao':mol.nao_nr(),'energy':mf.e_tot})
reference={'source':'PySCF '+pyscf.__version__,'guess':'1e','cases':cases}
(Path(__file__).parent/'fixtures/workspace.json').write_text(json.dumps(reference,indent=2)+'\n')
