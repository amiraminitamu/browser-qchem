"""Independent PySCF DF-JK comparisons; writes portable real-WASM fixtures."""
import ctypes as ct,json,sys,argparse
from pathlib import Path
import numpy as np
from pyscf import gto,scf,df,mp,lib
from check_molecular import load_library
from check_rhf import configure
lib.num_threads(1)
# Some managed runners do not expose process memory accounting through /proc.
try: lib.current_memory()
except FileNotFoundError: lib.current_memory=lambda:(0,0)
ROOT=Path(__file__).resolve().parents[1]
p=argparse.ArgumentParser();p.add_argument('library');p.add_argument('--write-fixtures',action='store_true');args=p.parse_args()
qc=configure(load_library(args.library));ip=ct.POINTER(ct.c_int);dp=ct.POINTER(ct.c_double);handle=ct.c_void_p
qc.qc_compute_mode.argtypes=[ct.c_int,*qc.qc_compute.argtypes];qc.qc_compute_mode.restype=handle
qc.qc_eri.argtypes=[handle,*([ct.c_int]*4)];qc.qc_eri.restype=ct.c_double
qc.qc_df_compute.argtypes=[ct.c_int,ip,dp,ct.c_int,ip,ct.c_int,dp,dp,ct.c_int,ip,ct.c_int,dp,dp];qc.qc_df_compute.restype=handle
for name in ['qc_df_error','qc_dfmp2_error','qc_uhf_error','qc_mp2_error']:getattr(qc,name).restype=ct.c_char_p
for name in ['qc_df_naux','qc_rhf_df_rank','qc_uhf_df_rank','qc_uhf_step','qc_uhf_status','qc_dfmp2_status','qc_dfmp2_step','qc_mp2_status','qc_mp2_step']:
 f=getattr(qc,name);f.argtypes=[handle];f.restype=ct.c_int
for name in ['qc_df_free','qc_uhf_free','qc_dfmp2_free','qc_mp2_free']:
 f=getattr(qc,name);f.argtypes=[handle];f.restype=None
for name in ['qc_df_data','qc_uhf_data']:
 f=getattr(qc,name);f.argtypes=[handle,ct.c_int];f.restype=dp
for name in ['qc_uhf_scalar','qc_dfmp2_scalar','qc_mp2_scalar']:
 f=getattr(qc,name);f.argtypes=[handle,ct.c_int];f.restype=ct.c_double
qc.qc_rhf_create_df.argtypes=[handle,handle,ct.c_int,ct.c_int,ct.c_int,*([ct.c_double]*3)];qc.qc_rhf_create_df.restype=handle
qc.qc_uhf_create_df.argtypes=[handle,handle,*([ct.c_int]*5),*([ct.c_double]*3)];qc.qc_uhf_create_df.restype=handle
for name in ['qc_dfmp2_rhf_create','qc_dfmp2_uhf_create']:
 f=getattr(qc,name);f.argtypes=[handle,handle];f.restype=handle
qc.qc_mp2_create.argtypes=[handle];qc.qc_mp2_create.restype=handle

def packed(mol):
 sh=[];ex=[];co=[]
 for a in range(mol.natm):
  for l,*ps in mol._basis[mol.atom_symbol(a)]:
   for k in range(1,len(ps[0])):sh.extend([a,l,len(ps)]);ex.extend(p[0] for p in ps);co.extend(p[k] for p in ps)
 return np.asarray(sh,dtype=np.int32),np.asarray(ex,dtype=np.float64),np.asarray(co,dtype=np.float64)
def inputs(mol):
 sh,ex,co=packed(mol);ch=np.ascontiguousarray(mol.atom_charges(),dtype=np.int32);xyz=np.ascontiguousarray(mol.atom_coords())
 return [mol.natm,ch.ctypes.data_as(ip),xyz.ctypes.data_as(dp),len(sh)//3,sh.ctypes.data_as(ip),len(ex),ex.ctypes.data_as(dp),co.ctypes.data_as(dp)],(ch,xyz,sh,ex,co)
def fitting(mol,basis):
 aux=df.addons.make_auxmol(mol,basis);values,keep=inputs(mol);sh,ex,co=packed(aux)
 h=qc.qc_df_compute(*values,len(sh)//3,sh.ctypes.data_as(ip),len(ex),ex.ctypes.data_as(dp),co.ctypes.data_as(dp));assert h,qc.qc_df_error();return h,aux
examples=json.loads((ROOT/'tests/fixtures/df-mp2.json').read_text())['cases'];fixtures=[]
for c in examples:
 mol=gto.M(atom=c['geometry'],unit=c['unit'],basis=c['basis'],charge=c['options']['charge'],spin=c['options']['multiplicity']-1,verbose=0)
 values,keep=inputs(mol);ints=qc.qc_compute_mode(0,*values);assert ints,qc.qc_error();assert not qc.qc_data(ints,3)
 fit,aux=fitting(mol,'def2-universal-jkfit');n=mol.nao_nr();m=aux.nao_nr()
 # Check every JKFIT 2c/3c entry independently, not just a fitted energy.
 metric=np.ctypeslib.as_array(qc.qc_df_data(fit,0),shape=(m*m,)).reshape(m,m)
 three=np.ctypeslib.as_array(qc.qc_df_data(fit,1),shape=(n*n*m,)).reshape(n,n,m)
 assert np.max(np.abs(metric-aux.intor('int2c2e')))<2e-9
 assert np.max(np.abs(three-df.incore.aux_e2(mol,aux,aosym='s1')))<2e-10
 u=c['options']['method']=='df-ump2';prefix='uhf' if u else 'rhf';ref=None;post=None;ri=None
 try:
  ref=qc.qc_uhf_create_df(ints,fit,*mol.nelec,500,1,int(c['options']['guess']=='mix'),1e-10,1e-8,1e-8) if u else qc.qc_rhf_create_df(ints,fit,mol.nelectron,500,1,1e-10,1e-8,1e-8)
  assert ref,getattr(qc,f'qc_{prefix}_error')()
  rank=getattr(qc,f'qc_{prefix}_df_rank')(ref)
  # Solver owns its whitened factors and survives early auxiliary handle release.
  qc.qc_df_free(fit);fit=None
  status=0
  while status==0:status=getattr(qc,f'qc_{prefix}_step')(ref)
  assert status==1,(c['name'],c['basis'],prefix,getattr(qc,f'qc_{prefix}_error')())
  read=lambda k,size:np.ctypeslib.as_array(getattr(qc,f'qc_{prefix}_data')(ref,k),shape=(size,)).copy()
  density=np.asarray([read(0,n*n).reshape(n,n),read(1,n*n).reshape(n,n)]) if u else read(0,n*n).reshape(n,n)
  mf=(scf.UHF(mol) if u else scf.RHF(mol)).density_fit(auxbasis='def2-universal-jkfit')
  mf.conv_tol=1e-12;mf.conv_tol_grad=1e-9;mf.kernel(dm0=density);assert mf.converged
  energy=getattr(qc,f'qc_{prefix}_scalar')(ref,0);assert abs(energy-mf.e_tot)<2e-9,(energy,mf.e_tot)
  fock=mf.get_fock(dm=density)
  native_f=np.asarray([read(2,n*n).reshape(n,n),read(3,n*n).reshape(n,n)]) if u else read(1,n*n).reshape(n,n)
  assert np.max(np.abs(native_f-fock))<2e-8
  s2=float(mf.spin_square()[0]) if u else 0
  if u:assert abs(qc.qc_uhf_scalar(ref,5)-s2)<2e-7
  eri=mol.intor('int2e');indices=[0,n-1,min(1,n-1),n-1]
  assert abs(qc.qc_eri(ints,*indices)-eri[tuple(indices)])<2e-10
  ri,riaux=fitting(mol,'cc-pvdz-ri');calc=mp.MP2(mf).density_fit(auxbasis='cc-pvdz-ri');calc.kernel(with_t2=False)
  post=(qc.qc_dfmp2_uhf_create if u else qc.qc_dfmp2_rhf_create)(ref,ri);assert post,qc.qc_dfmp2_error()
  while qc.qc_dfmp2_status(post)==0:qc.qc_dfmp2_step(post)
  assert qc.qc_dfmp2_status(post)==1,qc.qc_dfmp2_error()
  assert abs(qc.qc_dfmp2_scalar(post,1)-calc.e_corr)<2e-8
  qc.qc_dfmp2_free(post);post=None
  conv=None
  if not u:
   # Conventional MP2 on DF-SCF must use exact correlation ERIs without a dense AO tensor.
   conventional=mp.mp2.MP2(mf.undo_df())
   if n>mol.nelectron//2:conventional.kernel(with_t2=False)
   else:conventional.e_corr=0.0
   post=qc.qc_mp2_create(ref);assert post,qc.qc_mp2_error()
   while qc.qc_mp2_status(post)==0:qc.qc_mp2_step(post)
   assert qc.qc_mp2_status(post)==1,qc.qc_mp2_error()
   conv=float(conventional.e_corr);assert abs(qc.qc_mp2_scalar(post,1)-conv)<2e-8
   qc.qc_mp2_free(post);post=None
  assert not qc.qc_data(ints,3)
  fixtures.append({**c,'options':{**c['options'],'scfIntegrals':'df'},'referenceEnergy':float(mf.e_tot),'energy':float(calc.e_tot),'correlationEnergy':float(calc.e_corr),'oppositeSpinEnergy':float(calc.e_corr_os),'sameSpinEnergy':float(calc.e_corr_ss),'conventionalCorrelationEnergy':conv,'jkfitFunctions':m,'jkfitRank':rank,'spinSquare':s2,'eriIndices':indices,'exactERI':float(eri[tuple(indices)])})
  print(c['name'],c['basis'],prefix,'DF-JK PASS',flush=True)
 finally:
  if fit:qc.qc_df_free(fit)
  if ri:qc.qc_df_free(ri)
  if ref:getattr(qc,f'qc_{prefix}_free')(ref)
  qc.qc_free(ints)
if args.write_fixtures:(ROOT/'tests/fixtures/df-scf.json').write_text(json.dumps({'reference':'PySCF 2.14.0, def2-universal-JKFIT SCF + independent cc-pVDZ-RI MP2 or exact MP2','cases':fixtures},indent=2)+'\n')
print(len(fixtures),'DF-JK SCF and DF-MP2 comparisons, direct conventional MP2 and on-demand ERI PASS')
