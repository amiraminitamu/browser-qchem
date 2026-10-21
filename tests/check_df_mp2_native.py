"""Validate native 2c/3c integrals and both RI-MP2 spin variants against PySCF 2.14.
PYTHONPATH=/tmp/lorentzian-dev python tests/check_df_mp2_native.py build/native/libqc_core.so --write-fixtures
"""
import argparse,ctypes as ct,json
from pathlib import Path
import numpy as np
from pyscf import gto,scf,mp,df,lib
from check_molecular import load_library,native_handle
from check_rhf import configure
lib.num_threads(1)
ROOT=Path(__file__).resolve().parents[1]
parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('library');parser.add_argument('--write-fixtures',action='store_true');args=parser.parse_args()
qc=configure(load_library(args.library));ip=ct.POINTER(ct.c_int);dp=ct.POINTER(ct.c_double);handle=ct.c_void_p
qc.qc_df_compute.argtypes=[ct.c_int,ip,dp,ct.c_int,ip,ct.c_int,dp,dp,ct.c_int,ip,ct.c_int,dp,dp];qc.qc_df_compute.restype=handle
for name in ['qc_df_error','qc_dfmp2_error','qc_uhf_error']:getattr(qc,name).restype=ct.c_char_p
for name in ['qc_df_nao','qc_df_naux','qc_dfmp2_rank','qc_dfmp2_status','qc_dfmp2_step','qc_dfmp2_iteration','qc_dfmp2_steps','qc_uhf_step','qc_uhf_status']:
    f=getattr(qc,name);f.argtypes=[handle];f.restype=ct.c_int
for name in ['qc_df_free','qc_dfmp2_free','qc_uhf_free']:
    f=getattr(qc,name);f.argtypes=[handle];f.restype=None
for name in ['qc_df_data','qc_uhf_data']:
    f=getattr(qc,name);f.argtypes=[handle,ct.c_int];f.restype=dp
for name in ['qc_dfmp2_scalar','qc_uhf_scalar']:
    f=getattr(qc,name);f.argtypes=[handle,ct.c_int];f.restype=ct.c_double
for name in ['qc_dfmp2_rhf_create','qc_dfmp2_uhf_create']:
    f=getattr(qc,name);f.argtypes=[handle,handle];f.restype=handle
qc.qc_uhf_create.argtypes=[handle,*([ct.c_int]*5),*([ct.c_double]*3)];qc.qc_uhf_create.restype=handle

def packed(mol):
    shells=[];exps=[];coefs=[]
    for a in range(mol.natm):
        for l,*primitives in mol._basis[mol.atom_symbol(a)]:
            for k in range(1,len(primitives[0])):
                shells.extend([a,l,len(primitives)]);exps.extend(p[0] for p in primitives);coefs.extend(p[k] for p in primitives)
    return np.asarray(shells,dtype=np.int32),np.asarray(exps,dtype=np.float64),np.asarray(coefs,dtype=np.float64)
def df_handle(mol,aux):
    shell,exp,coef=packed(mol);ash,aexp,acoef=packed(aux)
    charges=np.ascontiguousarray(mol.atom_charges(),dtype=np.int32);xyz=np.ascontiguousarray(mol.atom_coords())
    result=qc.qc_df_compute(mol.natm,charges.ctypes.data_as(ip),xyz.ctypes.data_as(dp),len(shell)//3,shell.ctypes.data_as(ip),len(exp),exp.ctypes.data_as(dp),coef.ctypes.data_as(dp),len(ash)//3,ash.ctypes.data_as(ip),len(aexp),aexp.ctypes.data_as(dp),acoef.ctypes.data_as(dp))
    assert result,qc.qc_df_error().decode();return result

def compare_integrals(mol):
    aux=df.addons.make_auxmol(mol,'cc-pvdz-ri');fitting=df_handle(mol,aux)
    n,m=mol.nao_nr(),aux.nao_nr()
    metric=np.ctypeslib.as_array(qc.qc_df_data(fitting,0),shape=(m*m,)).reshape(m,m)
    three=np.ctypeslib.as_array(qc.qc_df_data(fitting,1),shape=(n*n*m,)).reshape(n,n,m)
    assert np.max(np.abs(metric-aux.intor('int2c2e')))<2e-9
    assert np.max(np.abs(three-df.incore.aux_e2(mol,aux,aosym='s1')))<2e-10
    return fitting,m

# Exercise normalization and angular ordering for all supported elements, including auxiliary f functions.
for charge,symbol in enumerate(['H','He','Li','Be','B','C','N','O','F','Ne'],1):
    atom=gto.M(atom=f'{symbol} .12 -.25 .31',basis='cc-pvdz',spin=charge%2,verbose=0)
    fitting,_=compare_integrals(atom);qc.qc_df_free(fitting)
print('All H-Ne orbital/auxiliary 2c and 3c integrals agree with PySCF PASS',flush=True)

examples=[('water','O .14 -.23 .31; H .79 .02 .91; H -.51 .61 .44',0,1),
          ('hydrogen','H 0 0 0; H 0 0 .74',0,1),('stretched_hydrogen','H 0 0 0; H 0 0 3',0,1),('hehplus','He 0 0 0; H 0 0 .8',1,1),
          ('methyl','C 0 0 0; H 1.07 0 0; H -.535 .92665 0; H -.535 -.92665 0',0,2),
          ('oxygen','O 0 0 0; O 0 0 1.21',0,3),('helium','He 0 0 0',0,1),
          ('h2plus','H 0 0 0; H 0 0 1.06',1,2),('hydrogen_atom','H 0 0 0',0,2)]
fixtures=[]
for basis in ['sto-3g','6-31g','cc-pvdz']:
  for name,geometry,charge,multiplicity in examples:
    for unrestricted in ([False,True] if multiplicity==1 else [True]):
      mol=gto.M(atom=geometry,basis=basis,charge=charge,spin=multiplicity-1,verbose=0)
      ints=native_handle(qc,mol);reference=post=fitting=None
      try:
        fitting,naux=compare_integrals(mol)
        guess='mix' if unrestricted and name=='stretched_hydrogen' else 'core'
        if unrestricted:reference=qc.qc_uhf_create(ints,*mol.nelec,500,1,int(guess=='mix'),1e-10,1e-8,1e-8)
        else:reference=qc.qc_rhf_create(ints,mol.nelectron,500,1,1e-10,1e-8,1e-8)
        prefix='uhf' if unrestricted else 'rhf';assert reference,getattr(qc,f'qc_{prefix}_error')().decode();status=0
        while status==0:status=getattr(qc,f'qc_{prefix}_step')(reference)
        assert status==1,(name,basis,prefix,getattr(qc,f'qc_{prefix}_error')().decode())
        n=mol.nao_nr();value=lambda kind:getattr(qc,f'qc_{prefix}_scalar')(reference,kind)
        data=lambda kind,size:np.ctypeslib.as_array(getattr(qc,f'qc_{prefix}_data')(reference,kind),shape=(size,)).copy()
        mf=scf.UHF(mol) if unrestricted else scf.RHF(mol)
        # Independent PySCF SCF, initialized at the native solution to choose the same stationary branch.
        density=np.asarray([data(0,n*n).reshape(n,n),data(1,n*n).reshape(n,n)]) if unrestricted else data(0,n*n).reshape(n,n)
        mf.conv_tol=1e-12;mf.conv_tol_grad=1e-9;mf.kernel(dm0=density);assert mf.converged
        assert abs(mf.e_tot-value(0))<2e-9,(name,basis,prefix,mf.e_tot,value(0))
        calc=mp.MP2(mf).density_fit(auxbasis='cc-pvdz-ri');calc.kernel(with_t2=False)
        post=(qc.qc_dfmp2_uhf_create if unrestricted else qc.qc_dfmp2_rhf_create)(reference,fitting)
        assert post,qc.qc_dfmp2_error().decode();status=qc.qc_dfmp2_status(post)
        while status==0:status=qc.qc_dfmp2_step(post)
        assert status==1,qc.qc_dfmp2_error().decode()
        expected={'energy':float(calc.e_tot),'correlationEnergy':float(calc.e_corr),'oppositeSpinEnergy':float(calc.e_corr_os),'sameSpinEnergy':float(calc.e_corr_ss)}
        for k,key in enumerate(expected):assert abs(qc.qc_dfmp2_scalar(post,k)-expected[key])<2e-8,(name,basis,prefix,key,qc.qc_dfmp2_scalar(post,k),expected[key])
        if mol.nelectron==1:assert qc.qc_dfmp2_scalar(post,1)==0
        assert qc.qc_dfmp2_rank(post)<=naux
        fixtures.append({'name':name,'geometry':geometry.replace(';','\n'),'unit':'angstrom','basis':basis,'options':{'method':'df-ump2' if unrestricted else 'df-rmp2','charge':charge,'multiplicity':multiplicity,'guess':guess},'referenceEnergy':float(mf.e_tot),'naux':naux,**expected})
        print(name,basis,'DF-UMP2' if unrestricted else 'DF-RMP2',f'{calc.e_corr:.12f}', 'PASS',flush=True)
      finally:
        if post:qc.qc_dfmp2_free(post)
        if fitting:qc.qc_df_free(fitting)
        if reference:getattr(qc,f'qc_{prefix}_free')(reference)
        qc.qc_free(ints)
if args.write_fixtures:(ROOT/'tests/fixtures/df-mp2.json').write_text(json.dumps({'reference':'PySCF 2.14.0, exact SCF + cc-pVDZ-RI correlation fitting','cases':fixtures},indent=2)+'\n')
print(f'{len(fixtures)} restricted/unrestricted DF-MP2 reference comparisons PASS')

# Redundant auxiliary functions exercise metric truncation; constructor owns all
# post-SCF state, so freeing its input handles before contraction must remain safe.
mol=gto.M(atom='O 0 0 0; H 0 -.757 .587; H 0 .757 .587',basis='sto-3g',verbose=0)
ints=native_handle(qc,mol);ref=qc.qc_rhf_create(ints,mol.nelectron,100,1,1e-10,1e-8,1e-8)
status=0
while status==0:status=qc.qc_rhf_step(ref)
assert status==1
base=df.addons.make_auxmol(mol,'cc-pvdz-ri');fitting=df_handle(mol,base)
post=qc.qc_dfmp2_rhf_create(ref,fitting);assert post
while qc.qc_dfmp2_status(post)==0:qc.qc_dfmp2_step(post)
baseline=qc.qc_dfmp2_scalar(post,1);qc.qc_dfmp2_free(post);qc.qc_df_free(fitting)
basis={s:gto.basis.load('cc-pvdz-ri',s)+[gto.basis.load('cc-pvdz-ri',s)[0]] for s in ['H','O']}
redundant=df.addons.make_auxmol(mol,basis);fitting=df_handle(mol,redundant)
metric=np.ctypeslib.as_array(qc.qc_df_data(fitting,0),shape=(redundant.nao_nr()**2,))
saved=metric[0];metric[0]=-1e6
assert not qc.qc_dfmp2_rhf_create(ref,fitting);assert b'positive semidefinite' in qc.qc_dfmp2_error()
metric[0]=saved
post=qc.qc_dfmp2_rhf_create(ref,fitting);assert post,qc.qc_dfmp2_error().decode()
assert qc.qc_dfmp2_rank(post)<redundant.nao_nr()
qc.qc_df_free(fitting);qc.qc_rhf_free(ref);qc.qc_free(ints)
while qc.qc_dfmp2_status(post)==0:qc.qc_dfmp2_step(post)
assert qc.qc_dfmp2_status(post)==1,qc.qc_dfmp2_error().decode()
assert abs(qc.qc_dfmp2_scalar(post,1)-baseline)<2e-9
qc.qc_dfmp2_free(post)
assert not qc.qc_dfmp2_rhf_create(None,None)
assert not qc.qc_dfmp2_uhf_create(None,None)
assert qc.qc_dfmp2_step(None)==-1
print('Redundant auxiliary metric, indefinite-metric rejection, owned-handle lifetime and recovery PASS')
