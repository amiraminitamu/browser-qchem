// SPDX-License-Identifier: AGPL-3.0-only
#include "df_mp2.h"
#include <Eigen/Dense>
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>
using Matrix=Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>;
using Vector=Eigen::VectorXd;
struct QCDFMP2 {
  int n,rank,o[2],v[2],iteration=0,status=0,steps;bool restricted;
  double reference,os=0,ss[2]={0,0},minimum=std::numeric_limits<double>::infinity(),gap,metric_min;
  std::array<Vector,2> eps;
  std::array<std::vector<Matrix>,2> factors;
};
static std::string error;
static constexpr double metric_cutoff=1e-10;
static std::unique_ptr<QCDFMP2> create(const QCDF *df,int n,bool restricted,double energy,
    const double *ca,const double *cb,const double *ea,const double *eb,int oa,int ob) {
  if(!df||qc_df_nao(df)!=n)throw std::runtime_error("DF and reference orbital dimensions differ.");
  auto r=std::make_unique<QCDFMP2>();r->n=n;r->restricted=restricted;r->reference=energy;
  r->o[0]=oa;r->o[1]=ob;r->v[0]=n-oa;r->v[1]=n-ob;r->steps=restricted?oa:oa+ob;
  r->eps[0]=Eigen::Map<const Vector>(ea,n);r->eps[1]=Eigen::Map<const Vector>(eb,n);
  r->gap=std::numeric_limits<double>::infinity();
  for(int s=0;s<2;++s)if(r->o[s]&&r->v[s])r->gap=std::min(r->gap,r->eps[s][r->o[s]]-r->eps[s][r->o[s]-1]);
  const int m=qc_df_naux(df);const Eigen::Map<const Matrix> metric(qc_df_data(df,0),m,m);
  Eigen::SelfAdjointEigenSolver<Matrix> eig((.5*(metric+metric.transpose())).eval());
  if(eig.info()!=Eigen::Success||!eig.eigenvalues().allFinite())throw std::runtime_error("DF Coulomb metric diagonalization failed.");
  r->metric_min=eig.eigenvalues()[0];
  if(r->metric_min < -std::max(1e-8,1e-10*eig.eigenvalues()[m-1]))throw std::runtime_error("DF Coulomb metric is not positive semidefinite.");
  r->rank=(eig.eigenvalues().array()>metric_cutoff).count();
  if(!r->rank)throw std::runtime_error("All auxiliary metric directions were removed.");
  Matrix whitening=eig.eigenvectors().rightCols(r->rank);
  for(int p=0;p<r->rank;++p)whitening.col(p)/=std::sqrt(eig.eigenvalues()[m-r->rank+p]);
  const Eigen::Map<const Matrix> three(qc_df_data(df,1),n*n,m);
  const Matrix b=three*whitening; // (mu nu|P) J^{-1/2}; storage O(N^2 M).
  if(!b.allFinite())throw std::runtime_error("DF whitening produced nonfinite factors.");
  for(int s=0;s<(restricted?1:2);++s){
    const Eigen::Map<const Matrix> c(s?cb:ca,n,n);
    for(int i=0;i<r->o[s];++i){
      Matrix x=Matrix::Zero(n,r->rank);
      for(int mu=0;mu<n;++mu)x.noalias()+=c(mu,i)*b.middleRows(mu*n,n);
      r->factors[s].push_back(c.rightCols(r->v[s]).transpose()*x);
    }
  }
  if((restricted&&!r->v[0])||(!restricted&&(oa+ob<2||(!r->v[0]&&!r->v[1]))))r->status=1;
  return r;
}
extern "C" QCDFMP2 *qc_dfmp2_rhf_create(const QCRHF *ref,const QCDF *df){
  error.clear();try{
    if(!ref||qc_rhf_status(ref)!=1)throw std::runtime_error("DF-RMP2 requires a converged RHF reference.");
    const int n=qc_nao(qc_rhf_integrals(ref));int o=0;for(int i=0;i<n;++i)if(qc_rhf_data(ref,4)[i]>1)++o;
    return create(df,n,true,qc_rhf_scalar(ref,0),qc_rhf_data(ref,2),qc_rhf_data(ref,2),qc_rhf_data(ref,3),qc_rhf_data(ref,3),o,o).release();
  }catch(const std::exception &e){error=e.what();return nullptr;}
}
extern "C" QCDFMP2 *qc_dfmp2_uhf_create(const QCUHF *ref,const QCDF *df){
  error.clear();try{
    if(!ref||qc_uhf_status(ref)!=1)throw std::runtime_error("DF-UMP2 requires a converged UHF reference.");
    const int n=qc_uhf_nao(ref);int oa=0,ob=0;
    for(int i=0;i<n;++i){oa+=qc_uhf_data(ref,8)[i]>.5;ob+=qc_uhf_data(ref,9)[i]>.5;}
    return create(df,n,false,qc_uhf_scalar(ref,0),qc_uhf_data(ref,4),qc_uhf_data(ref,5),qc_uhf_data(ref,6),qc_uhf_data(ref,7),oa,ob).release();
  }catch(const std::exception &e){error=e.what();return nullptr;}
}
static double denominator(QCDFMP2 &r,int s,int t,int i,int j,int a,int b){
  const double d=r.eps[s][i]+r.eps[t][j]-r.eps[s][r.o[s]+a]-r.eps[t][r.o[t]+b];
  if(!std::isfinite(d)||d>=-1e-10)throw std::runtime_error("DF-MP2 denominator is non-negative or too small. The reference is unsuitable for MP2.");
  r.minimum=std::min(r.minimum,-d);return d;
}
extern "C" int qc_dfmp2_step(QCDFMP2 *r){
  if(!r){error="Invalid DF-MP2 handle.";return -1;}if(r->status)return r->status;
  try{
    const int s=r->restricted||r->iteration<r->o[0]?0:1,i=r->iteration-(s?r->o[0]:0),v=r->v[s];
    double os=0,ss=0;
    for(int j=0;j<r->o[s];++j){
      if(!r->restricted&&i==j)continue; // No same-spin excitation from the same occupied orbital.
      const Matrix g=r->factors[s][i]*r->factors[s][j].transpose();
      for(int a=0;a<v;++a)for(int b=0;b<v;++b){
        if(!r->restricted&&a==b)continue;
        const double d=denominator(*r,s,s,i,j,a,b),direct=g(a,b),exchange=g(b,a);
        if(r->restricted)os+=direct*direct/d;
        ss+=(r->restricted?1:.5)*direct*(direct-exchange)/d;
      }
    }
    if(!r->restricted&&s==0)for(int j=0;j<r->o[1];++j){
      const Matrix g=r->factors[0][i]*r->factors[1][j].transpose();
      for(int a=0;a<r->v[0];++a)for(int b=0;b<r->v[1];++b)os+=g(a,b)*g(a,b)/denominator(*r,0,1,i,j,a,b);
    }
    if(!std::isfinite(os+ss))throw std::runtime_error("DF-MP2 correlation energy is not finite.");
    r->os+=os;
    if(r->restricted){r->ss[0]+=.5*ss;r->ss[1]+=.5*ss;}else r->ss[s]+=ss;
    if(++r->iteration==r->steps)r->status=1;
  }catch(const std::exception &e){error=e.what();r->status=-1;}
  return r->status;
}
extern "C" const char *qc_dfmp2_error(){return error.c_str();}
extern "C" int qc_dfmp2_status(const QCDFMP2 *r){return r?r->status:-1;}
extern "C" int qc_dfmp2_iteration(const QCDFMP2 *r){return r?r->iteration:0;}
extern "C" int qc_dfmp2_steps(const QCDFMP2 *r){return r?r->steps:0;}
extern "C" int qc_dfmp2_nocc(const QCDFMP2 *r,int spin){return r&&spin>=0&&spin<2?r->o[spin]:0;}
extern "C" int qc_dfmp2_nvirt(const QCDFMP2 *r,int spin){return r&&spin>=0&&spin<2?r->v[spin]:0;}
extern "C" int qc_dfmp2_rank(const QCDFMP2 *r){return r?r->rank:0;}
extern "C" double qc_dfmp2_scalar(const QCDFMP2 *r,int kind){
  if(r)switch(kind){case 0:return r->reference+r->os+r->ss[0]+r->ss[1];case 1:return r->os+r->ss[0]+r->ss[1];case 2:return r->os;case 3:return r->ss[0]+r->ss[1];case 4:return r->minimum;case 5:return r->gap;case 6:return r->reference;case 7:return r->ss[0];case 8:return r->ss[1];case 9:return metric_cutoff;case 10:return r->metric_min;}
  return std::numeric_limits<double>::quiet_NaN();
}
extern "C" void qc_dfmp2_free(QCDFMP2 *r){delete r;}
