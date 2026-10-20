"""Compare real WASM UHF and conventional all-electron RHF-MP2 with PySCF.
Run from the repository: python3 tests/check_methods.py
"""
import json,subprocess
from pathlib import Path
import numpy as np
from pyscf import gto,scf,mp,lib
lib.num_threads(1)
ROOT=Path(__file__).resolve().parents[1]
WATER='O 0 0 0; H 0 -.757 .587; H 0 .757 .587'
CASES=[]
for basis in ['sto-3g','6-31g','cc-pvdz']:
    for name,geom,charge,spin,guess in [('H','H 0 0 0',0,1,'core'),('He+','He 0 0 0',1,1,'core'),('H2+','H 0 0 0; H 0 0 1.06',1,1,'core'),('HeH','He 0 0 0; H 0 0 1.0',0,1,'core'),('HeH+','He 0 0 0; H 0 0 .8',1,0,'core'),('H-','H 0 0 0',-1,0,'core'),('Li','Li 0 0 0',0,1,'core'),('OH','O 0 0 0; H .1 .2 .94',0,1,'core'),('water',WATER,0,0,'core'),('water mixed',WATER,0,0,'mix'),('H2 equilibrium','H 0 0 0; H 0 0 .74',0,0,'mix'),('H2 broken','H 0 0 0; H 0 0 3',0,0,'mix')]:
        CASES.append((name,geom,basis,charge,spin,guess,'uhf'))
    for name,geom in [('water',WATER),('HeH+','He 0 0 0; H 0 0 .8'),('H2','H 0 0 0; H 0 0 .74'),('HF','F .1 .2 .3; H .3 .1 1.2'),('He','He 0 0 0'),('ammonia','N 0 0 .1; H .93 0 -.25; H -.465 .8054 -.25; H -.465 -.8054 -.25')]:
        CASES.append((name,geom,basis,1 if name=='HeH+' else 0,0,'core','mp2'))
CASES.extend([('O triplet','O 0 0 0','sto-3g',0,2,'core','uhf'),('O2 triplet','O 0 0 0; O 0 0 1.21','sto-3g',0,2,'core','uhf')])
references=[]
for name,geom,basis,charge,spin,guess,method in CASES:
    mol=gto.M(atom=geom,basis=basis,charge=charge,spin=spin,verbose=0);mol.incore_anyway=True
    mf=scf.UHF(mol) if method=='uhf' else scf.RHF(mol)
    mf.conv_tol=1e-12;mf.conv_tol_grad=1e-10;mf.max_cycle=500
    mf.init_guess='1e';mf.init_guess_breaksym=0
    request={'geometry':geom.replace(';','\n'),'unit':'angstrom','basis':basis,'options':{'method':method,'charge':charge,'multiplicity':spin+1,'guess':guess,'maxIterations':500}}
    if method=='uhf':
        # Match the exact starting density to avoid different stationary roots
        # caused by degenerate core orbitals and eigensolver-dependent rotations.
        initial=json.loads(subprocess.run(['node','tests/export_methods.mjs','--initial-uhf'],cwd=ROOT,input=json.dumps(request),text=True,check=True,capture_output=True).stdout)
        mf.kernel(dm0=np.array(initial).reshape(2,mol.nao_nr(),mol.nao_nr()))
    else:mf.kernel()
    assert mf.converged,(name,basis)
    result=json.loads(subprocess.run(['node','tests/export_methods.mjs'],cwd=ROOT,input=json.dumps(request),text=True,check=True,capture_output=True).stdout)
    r=result['uhf' if method=='uhf' else 'rhf'];assert r['converged'],(name,basis,r['history'][-1])
    initial_root_energy=float(mf.e_tot)
    matched_stationary_root=False
    if method=='uhf' and abs(r['energy']-mf.e_tot)>=2e-8:
        # DIIS trajectories can select distinct stationary states even from the
        # same core guess. Independently reconverge and verify the reported root.
        # This does not claim that either trajectory finds the ground state.
        dm=np.array([r['densityAlpha'],r['densityBeta']]).reshape(2,mol.nao_nr(),mol.nao_nr())
        mf.kernel(dm0=dm);assert mf.converged
        matched_stationary_root=True
        print(f'  {name}/{basis}: distinct SCF root; independently reconverged reported density with PySCF',flush=True)
    error=abs(r['energy']-mf.e_tot);assert error<2e-8,(name,basis,error,r['energy'],mf.e_tot)
    ref={**request,'name':name,'referenceEnergy':float(mf.e_tot),'matchedStationaryRoot':matched_stationary_root,'pyscfInitialTrajectoryEnergy':initial_root_energy}
    if method=='uhf':
        assert abs(r['spinSquare']-mf.spin_square()[0])<2e-6,(name,basis,r['spinSquare'],mf.spin_square())
        S=mol.intor('int1e_ovlp');H=mf.get_hcore();eri=mol.intor('int2e');ds=[]
        for label,nocc in [('Alpha',mol.nelec[0]),('Beta',mol.nelec[1])]:
            n=mol.nao_nr();D=np.array(r['density'+label]).reshape(n,n);C=np.array(r['coefficients'+label]).reshape(n,n);F=np.array(r['fock'+label]).reshape(n,n);eps=np.array(r['orbitalEnergies'+label]);ds.append(D)
            assert np.max(abs(C.T@S@C-np.eye(n)))<2e-8
            assert abs(np.trace(D@S)-nocc)<1e-8
            assert np.max(abs(D@S@D-D))<2e-8
            assert np.max(abs(F@C-S@C@np.diag(eps)))<2e-6
        for label,D in zip(['Alpha','Beta'],ds):
            F=H+np.einsum('kl,ijkl->ij',ds[0]+ds[1],eri)-np.einsum('kl,ikjl->ij',D,eri)
            assert np.max(abs(F-np.array(r['fock'+label]).reshape(n,n)))<2e-8
        ref['spinSquare']=mf.spin_square()[0]
    else:
        pt=mp.MP2(mf);pt.kernel();m=result['mp2']
        for key,expected in [('energy',pt.e_tot),('correlationEnergy',pt.e_corr),('oppositeSpinEnergy',pt.e_corr_os),('sameSpinEnergy',pt.e_corr_ss)]:
            assert abs(m[key]-expected)<2e-8,(name,basis,key,m[key],expected)
        ref.update(energy=pt.e_tot,correlationEnergy=pt.e_corr,oppositeSpinEnergy=pt.e_corr_os,sameSpinEnergy=pt.e_corr_ss)
    references.append(ref)
    print(f'{method.upper()} {name}/{basis}: energy, '+('spin and SCF identities' if method=='uhf' else 'OS/SS correlation')+' PASS',flush=True)
ion_scans=[]
for basis in ['sto-3g','6-31g','cc-pvdz']:
    for name,first,charge,multiplicity in [('H2+','H',1,2),('HeH','He',0,2),('HeH+','He',1,1)]:
        request={'name':name,'geometry':f'{first} 0 0 0\nH 0 0 2','unit':'bohr','basis':basis,'options':{'method':'uhf','charge':charge,'multiplicity':multiplicity,'guess':'core','maxIterations':500},'scan':{'a':0,'b':1,'moving':[1],'start':1.4,'end':4,'count':3,'unit':'bohr'}}
        actual=json.loads(subprocess.run(['node','tests/export_scan.mjs'],cwd=ROOT,input=json.dumps(request),text=True,check=True,capture_output=True).stdout)
        energies=[]
        for point in actual['points']:
            assert point['converged'],point
            mol=gto.M(atom=f'{first} 0 0 0; H 0 0 {point["coordinate"]}',unit='bohr',basis=basis,charge=charge,spin=multiplicity-1,verbose=0);mol.incore_anyway=True
            mf=scf.UHF(mol);mf.init_guess='1e';mf.init_guess_breaksym=0;mf.conv_tol=1e-12;mf.conv_tol_grad=1e-10;mf.max_cycle=500;mf.kernel();assert mf.converged
            assert abs(mf.e_tot-point['energy'])<2e-8,(name,basis,point['coordinate'],mf.e_tot,point['energy'])
            assert point['result']['uhf']['electrons']==mol.nelectron
            energies.append(float(mf.e_tot))
        ion_scans.append({**request,'energies':energies})
        print(f'UHF {name}/{basis}: charged/open-shell three-point PES PASS',flush=True)
(ROOT/'tests/fixtures/methods.json').write_text(json.dumps({'source':'PySCF '+__import__('pyscf').__version__,'cases':references,'ionScans':ion_scans},indent=2)+'\n')
print(f'{len(references)} single points and {len(ion_scans)*3} ion/open-shell scan points PASS')
