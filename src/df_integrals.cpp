// SPDX-License-Identifier: AGPL-3.0-only
#include "df_integrals.h"
#include "limits.h"
extern "C" {
#include "cint.h"
#include "cint_funcs.h"
// These built-in Coulomb operators are not listed in the generated 1e/2e header.
extern CINTIntegralFunction int2c2e_sph;
extern CINTIntegralFunction int3c2e_sph;
}
#include <algorithm>
#include <cmath>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>
struct QCDF {int n,m;std::vector<double> metric,three;};
static std::string error;
static void require(bool ok,const char *message){if(!ok)throw std::runtime_error(message);}
extern "C" QCDF *qc_df_compute(int natm,const int *charges,const double *xyz,
    int nshell,const int *shells,int nprim,const double *exponents,const double *coefficients,
    int nauxshell,const int *auxshells,int nauxprim,const double *auxexponents,const double *auxcoefficients) {
  error.clear();
  try {
    require(natm>=1&&natm<=32&&charges&&xyz,"Invalid DF atom arrays.");
    require(nshell>=1&&nshell<=QC_MAX_SHELLS&&shells&&nprim>=1&&nprim<=4096&&exponents&&coefficients,"Invalid DF orbital basis.");
    require(nauxshell>=1&&nauxshell<=1024&&auxshells&&nauxprim>=1&&nauxprim<=32768&&auxexponents&&auxcoefficients,"Invalid DF auxiliary basis.");
    const int count=nshell+nauxshell;
    std::vector<int> atoms(natm*ATM_SLOTS,0),basis(count*BAS_SLOTS,0),offsets(count+1,0);
    std::vector<double> env(PTR_ENV_START+3*natm+2*(nprim+nauxprim),0);
    for(int a=0;a<natm;++a){
      require(charges[a]>=1&&charges[a]<=118,"Invalid DF nuclear charge.");
      atoms[a*ATM_SLOTS+CHARGE_OF]=charges[a];atoms[a*ATM_SLOTS+PTR_COORD]=PTR_ENV_START+3*a;atoms[a*ATM_SLOTS+NUC_MOD_OF]=POINT_NUC;
      for(int d=0;d<3;++d){require(std::isfinite(xyz[3*a+d])&&std::abs(xyz[3*a+d])<=1e6,"Invalid DF coordinate.");env[PTR_ENV_START+3*a+d]=xyz[3*a+d];}
    }
    int position=PTR_ENV_START+3*natm,maxdim=1;
    for(int group=0;group<2;++group){
      const int *records=group?auxshells:shells;const double *exp=group?auxexponents:exponents,*coef=group?auxcoefficients:coefficients;
      const int first=group?nshell:0,number=group?nauxshell:nshell,primitives=group?nauxprim:nprim;int used=0;
      for(int s=0;s<number;++s){
        const int atom=records[3*s],l=records[3*s+1],np=records[3*s+2],k=first+s;
        require(atom>=0&&atom<natm&&l>=0&&l<=(group?6:2)&&np>=1&&np<=64&&used+np<=primitives,"Invalid DF shell record.");
        offsets[k+1]=offsets[k]+2*l+1;maxdim=std::max(maxdim,2*l+1);
        basis[k*BAS_SLOTS+ATOM_OF]=atom;basis[k*BAS_SLOTS+ANG_OF]=l;basis[k*BAS_SLOTS+NPRIM_OF]=np;basis[k*BAS_SLOTS+NCTR_OF]=1;
        basis[k*BAS_SLOTS+PTR_EXP]=position;basis[k*BAS_SLOTS+PTR_COEFF]=position+np;
        for(int p=0;p<np;++p)require(std::isfinite(exp[used+p])&&exp[used+p]>0&&exp[used+p]<=1e8&&std::isfinite(coef[used+p])&&std::abs(coef[used+p])<=1e8,"Invalid DF primitive.");
        double norm=0;
        for(int p=0;p<np;++p)for(int q=0;q<np;++q)norm+=coef[used+p]*coef[used+q]*std::pow(2*std::sqrt(exp[used+p]*exp[used+q])/(exp[used+p]+exp[used+q]),l+1.5);
        require(std::isfinite(norm)&&norm>1e-30,"Invalid DF contraction normalization.");
        for(int p=0;p<np;++p){env[position+p]=exp[used+p];env[position+np+p]=coef[used+p]*CINTgto_norm(l,exp[used+p])/std::sqrt(norm);}
        position+=2*np;used+=np;
      }
      require(used==primitives,"DF primitive count mismatch.");
    }
    auto r=std::make_unique<QCDF>();r->n=offsets[nshell];r->m=offsets[count]-r->n;
    require(r->n<=QC_MAX_AO&&r->m<=QC_MAX_AUXILIARY,"DF size exceeds the 256 orbital / 2048 auxiliary function limits.");
    const size_t n=r->n,m=r->m;
    require(8.0*(2.0*n*n*m+6.0*m*m+48.0*n*n)<=QC_MAX_FITTING_BYTES,"DF fitting workspace exceeds 384 MiB. Use a smaller basis or molecule.");
    r->metric.resize(m*m);r->three.resize(n*n*m);
    std::vector<double> buffer(static_cast<size_t>(maxdim)*maxdim*maxdim);
    // libcint's shell blocks are first-index-fastest; all public arrays are row-major.
    for(int p=nshell;p<count;++p)for(int q=nshell;q<count;++q){
      int shls[2]={p,q};const int dp=offsets[p+1]-offsets[p],dq=offsets[q+1]-offsets[q];
      std::fill(buffer.begin(),buffer.end(),0);int2c2e_sph(buffer.data(),nullptr,shls,atoms.data(),natm,basis.data(),count,env.data(),nullptr,nullptr);
      for(int j=0;j<dq;++j)for(int i=0;i<dp;++i)r->metric[(offsets[p]-n+i)*m+offsets[q]-n+j]=buffer[i+dp*j];
    }
    for(int a=0;a<nshell;++a)for(int b=0;b<nshell;++b)for(int p=nshell;p<count;++p){
      int shls[3]={a,b,p};const int da=offsets[a+1]-offsets[a],db=offsets[b+1]-offsets[b],dp=offsets[p+1]-offsets[p];
      std::fill(buffer.begin(),buffer.end(),0);int3c2e_sph(buffer.data(),nullptr,shls,atoms.data(),natm,basis.data(),count,env.data(),nullptr,nullptr);
      for(int k=0;k<dp;++k)for(int j=0;j<db;++j)for(int i=0;i<da;++i)r->three[((offsets[a]+i)*n+offsets[b]+j)*m+offsets[p]-n+k]=buffer[i+da*(j+db*k)];
    }
    require(std::all_of(r->metric.begin(),r->metric.end(),[](double x){return std::isfinite(x);})&&std::all_of(r->three.begin(),r->three.end(),[](double x){return std::isfinite(x);}),"Nonfinite DF integral.");
    return r.release();
  }catch(const std::exception &e){error=e.what();return nullptr;}
}
extern "C" const char *qc_df_error(){return error.c_str();}
extern "C" int qc_df_nao(const QCDF *r){return r?r->n:0;}
extern "C" int qc_df_naux(const QCDF *r){return r?r->m:0;}
extern "C" const double *qc_df_data(const QCDF *r,int kind){return !r?nullptr:kind==0?r->metric.data():kind==1?r->three.data():nullptr;}
extern "C" void qc_df_free(QCDF *r){delete r;}
