"""Compare all integral entries and RHF results for all three bundled bases."""
import json,subprocess
from pathlib import Path
import numpy as np
from pyscf import gto,scf,lib
lib.num_threads(1)
ROOT=Path(__file__).resolve().parents[1]
cases=[('water','O 0 0 0; H 0 -.757 .587; H 0 .757 .587'),
       ('skew water','O .14 -.23 .31; H .79 .02 .91; H -.51 .61 .44'),
       ('ammonia','N 0 0 .1; H .93 0 -.25; H -.465 .8054 -.25; H -.465 -.8054 -.25'),
       ('methane','C 0 0 0; H .63 .63 .63; H -.63 -.63 .63; H -.63 .63 -.63; H .63 -.63 -.63')]
largest_integral=largest_energy=largest_density=0.0
for basis in ['sto-3g','6-31g','cc-pvdz']:
    for name,geometry in cases:
        mol=gto.M(atom=geometry,basis=basis,unit='Angstrom',verbose=0)
        mol.incore_anyway=True
        process=subprocess.run(['node',str(ROOT/'tests/export_rhf.mjs')],
          input=json.dumps({'geometry':geometry.replace(';','\n'),'unit':'angstrom','basis':basis}),
          text=True,capture_output=True,check=True,cwd=ROOT)
        actual=json.loads(process.stdout)
        assert actual['nao']==mol.nao_nr()
        errors=[]
        for key,op in [('overlap','int1e_ovlp_sph'),('kinetic','int1e_kin_sph'),('nuclear','int1e_nuc_sph'),('eri','int2e_sph')]:
            expected=mol.intor(op,aosym='s1').reshape(-1)
            error=float(np.max(np.abs(np.asarray(actual[key])-expected)))
            assert error<1e-10,(basis,name,key,error)
            errors.append(error)
        ref=scf.RHF(mol); ref.conv_tol=1e-12; ref.conv_tol_grad=1e-10; ref.init_guess='1e'; ref.max_cycle=200
        ref.kernel(); assert ref.converged
        assert actual['rhf']['converged'],(basis,name,'SCF failed')
        energy_error=abs(actual['rhf']['energy']-ref.e_tot)
        density=np.asarray(actual['rhf']['density']).reshape(mol.nao_nr(),mol.nao_nr())
        density_error=float(np.max(np.abs(density-ref.make_rdm1())))
        assert energy_error<1e-9,(basis,name,'energy',energy_error)
        assert density_error<1e-7,(basis,name,'density',density_error)
        largest_integral=max(largest_integral,*errors); largest_energy=max(largest_energy,energy_error); largest_density=max(largest_density,density_error)
        print(f'{basis:8s} {name:10s} {mol.nao_nr():2d} AOs | integral {max(errors):.2e}, energy {energy_error:.2e}, density {density_error:.2e} PASS',flush=True)
print(f'Maximum errors: integrals {largest_integral:.3e}, energy {largest_energy:.3e}, density {largest_density:.3e}')
for basis in ['sto-3g','6-31g','cc-pvdz']:
    maximum=0.0
    for charge,symbol in enumerate(['H','He','Li','Be','B','C','N','O','F','Ne'],1):
        mol=gto.M(atom=f'{symbol} .13 -.21 .37',basis=basis,spin=charge%2,verbose=0)
        mol.incore_anyway=True
        process=subprocess.run(['node',str(ROOT/'tests/export_integrals.mjs')],
          input=json.dumps({'geometry':f'{symbol} .13 -.21 .37','basis':basis,'unit':'angstrom'}),
          text=True,capture_output=True,check=True,cwd=ROOT)
        actual=json.loads(process.stdout)
        assert actual['nao']==mol.nao_nr()
        for key,op in [('overlap','int1e_ovlp_sph'),('kinetic','int1e_kin_sph'),('nuclear','int1e_nuc_sph'),('eri','int2e_sph')]:
            error=float(np.max(np.abs(np.asarray(actual[key])-mol.intor(op,aosym='s1').reshape(-1))))
            assert error<1e-10,(basis,symbol,key,error)
            maximum=max(maximum,error)
    print(f'{basis}: all H-Ne atomic integral entries, max error {maximum:.3e} PASS',flush=True)
