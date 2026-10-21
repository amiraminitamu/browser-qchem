// SPDX-License-Identifier: AGPL-3.0-only
#include "uhf.h"
#include "df_factors.h"
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
using Pair=std::array<Matrix,2>;
static std::string last_error;
struct QCUHF {
    int n,occ[2],maxiter,iteration=0,status=0;
    bool diis,mixed=false;
    double energy,enuc,de=0,drms=0,residual=0,overlap_min,etol,dtol,rtol;
    std::unique_ptr<DFFactors> df;
    const double *eri;
    Matrix s,h,x;
    Pair d,f,c;
    std::array<Vector,2> eps,occupations;
    std::vector<Pair> f_history,e_history;
};
static Pair fock(const QCUHF &r,const Pair &d) {
    Pair f={r.h,r.h};const size_t n=r.n;
    if(r.df){const Matrix j=r.df->coulomb((d[0]+d[1]).eval());for(int s=0;s<2;++s)f[s]+=j-r.df->exchange(r.c[s].leftCols(r.occ[s]));return f;}
    for(size_t i=0;i<n;++i)for(size_t j=0;j<n;++j){
        double coulomb=0,exchange[2]={0,0};
        for(size_t k=0;k<n;++k)for(size_t l=0;l<n;++l){
            coulomb+=(d[0](k,l)+d[1](k,l))*r.eri[((i*n+j)*n+k)*n+l];
            const double value=r.eri[((i*n+k)*n+j)*n+l];
            exchange[0]+=d[0](k,l)*value;exchange[1]+=d[1](k,l)*value;
        }
        for(int spin=0;spin<2;++spin)f[spin](i,j)+=coulomb-exchange[spin];
    }
    for(auto &matrix:f)matrix=(.5*(matrix+matrix.transpose())).eval();
    return f;
}
static double energy(const QCUHF &r,const Pair &d,const Pair &f){
    return r.enuc+.5*((d[0].array()*(r.h+f[0]).array()).sum()+(d[1].array()*(r.h+f[1]).array()).sum());
}
static Pair errors(const QCUHF &r,const Pair &d,const Pair &f){
    Pair e;for(int spin=0;spin<2;++spin)e[spin]=r.x.transpose()*(f[spin]*d[spin]*r.s-r.s*d[spin]*f[spin])*r.x;
    return e;
}
static double norm(const Pair &p){return std::sqrt(p[0].squaredNorm()+p[1].squaredNorm());}
static double difference(const Pair &a,const Pair &b){return std::sqrt((a[0]-b[0]).squaredNorm()+(a[1]-b[1]).squaredNorm());}
static Pair diagonalize(QCUHF &r,const Pair &f){
    Pair d;
    for(int spin=0;spin<2;++spin){
        Matrix transformed=r.x.transpose()*f[spin]*r.x;
        transformed=(.5*(transformed+transformed.transpose())).eval();
        Eigen::SelfAdjointEigenSolver<Matrix> eig(transformed);
        if(eig.info()!=Eigen::Success)throw std::runtime_error("UHF Fock diagonalization failed.");
        r.c[spin]=r.x*eig.eigenvectors();r.eps[spin]=eig.eigenvalues();
        const Matrix occupied=r.c[spin].leftCols(r.occ[spin]);d[spin]=occupied*occupied.transpose();
    }
    return d;
}
static Pair extrapolate(QCUHF &r){
    if(!r.diis)return r.f;
    r.f_history.push_back(r.f);r.e_history.push_back(errors(r,r.d,r.f));
    if(r.f_history.size()>8){r.f_history.erase(r.f_history.begin());r.e_history.erase(r.e_history.begin());}
    const int m=static_cast<int>(r.f_history.size());if(m<2)return r.f;
    Matrix b=Matrix::Zero(m+1,m+1);
    for(int i=0;i<m;++i)for(int j=0;j<m;++j)
        b(i,j)=(r.e_history[i][0].array()*r.e_history[j][0].array()).sum()+(r.e_history[i][1].array()*r.e_history[j][1].array()).sum();
    const double scale=b.topLeftCorner(m,m).cwiseAbs().maxCoeff();if(scale<1e-28)return r.f;
    b.topLeftCorner(m,m)/=scale;b.block(0,m,m,1).setConstant(-1);b.block(m,0,1,m).setConstant(-1);
    Vector rhs=Vector::Zero(m+1);rhs(m)=-1;
    Eigen::FullPivLU<Matrix> solve(b);solve.setThreshold(1e-12);if(!solve.isInvertible())return r.f;
    Vector weights=solve.solve(rhs);
    if(!weights.allFinite()||weights.head(m).cwiseAbs().maxCoeff()>1e6||(b*weights-rhs).norm()>1e-7)return r.f;
    Pair f={Matrix::Zero(r.n,r.n),Matrix::Zero(r.n,r.n)};
    for(int i=0;i<m;++i)for(int spin=0;spin<2;++spin)f[spin]+=weights(i)*r.f_history[i][spin];
    return f;
}
extern "C" const char *qc_uhf_error(void){return last_error.c_str();}
extern "C" QCUHF *qc_uhf_create_df(const QCResult *integrals,const QCDF *df,int na,int nb,int maxiter,int diis,int guess,double etol,double dtol,double rtol){
    last_error.clear();
    try{
        const int n=qc_nao(integrals);
        if(n<1||n>40||na<nb||nb<0||na<1||na>n||nb>n)throw std::runtime_error("UHF spin populations must fit the basis, with Nalpha >= Nbeta >= 0 and at least one electron.");
        if(maxiter<1||maxiter>500||guess<0||guess>1||!std::isfinite(etol)||!std::isfinite(dtol)||!std::isfinite(rtol)||etol<=0||dtol<=0||rtol<=0)
            throw std::runtime_error("Invalid UHF guess, iteration limit or convergence threshold.");
        auto r=std::make_unique<QCUHF>();r->n=n;r->occ[0]=na;r->occ[1]=nb;r->maxiter=maxiter;r->diis=diis!=0;
        r->etol=etol;r->dtol=dtol;r->rtol=rtol;r->enuc=qc_enuc(integrals);r->eri=qc_data(integrals,3);
        if(df)r->df=std::make_unique<DFFactors>(df,n);else if(!r->eri)throw std::runtime_error("Exact SCF requires four-center integrals.");
        r->s=Eigen::Map<const Matrix>(qc_data(integrals,0),n,n);
        r->h=Eigen::Map<const Matrix>(qc_data(integrals,1),n,n)+Eigen::Map<const Matrix>(qc_data(integrals,2),n,n);
        Eigen::SelfAdjointEigenSolver<Matrix> eig(r->s);
        if(eig.info()!=Eigen::Success)throw std::runtime_error("Overlap diagonalization failed.");
        r->overlap_min=eig.eigenvalues().minCoeff();
        if(r->overlap_min<1e-8)throw std::runtime_error("Overlap is nearly linearly dependent (smallest eigenvalue < 1e-8). Adjust the geometry.");
        r->x=eig.eigenvectors()*eig.eigenvalues().cwiseSqrt().cwiseInverse().asDiagonal()*eig.eigenvectors().transpose();
        r->d=diagonalize(*r,Pair{r->h,r->h});
        if(guess==1&&na==nb&&na>0&&na<n){
            for(int spin=0;spin<2;++spin){
                const Vector homo=r->c[spin].col(na-1),lumo=r->c[spin].col(na);
                const double sign=spin==0?1.0:-1.0;
                r->c[spin].col(na-1)=(homo+sign*lumo)/std::sqrt(2.0);
                r->c[spin].col(na)=(lumo-sign*homo)/std::sqrt(2.0);
                const Matrix occupied=r->c[spin].leftCols(na);r->d[spin]=occupied*occupied.transpose();
            }
            r->mixed=true;
        }
        for(int spin=0;spin<2;++spin){r->occupations[spin]=Vector::Zero(n);r->occupations[spin].head(r->occ[spin]).setConstant(1);}
        r->f=fock(*r,r->d);r->energy=energy(*r,r->d,r->f);r->residual=norm(errors(*r,r->d,r->f))/(n*std::sqrt(2.0));
        return r.release();
    }catch(const std::exception &e){last_error=e.what();return nullptr;}
}
extern "C" int qc_uhf_step(QCUHF *r){
    if(!r){last_error="Missing UHF handle.";return -1;}if(r->status!=0)return r->status;
    try{
        const double scale=r->n*std::sqrt(2.0);
        Pair next=diagonalize(*r,extrapolate(*r)),f=fock(*r,next);
        double e=energy(*r,next,f);r->de=e-r->energy;r->drms=difference(next,r->d)/scale;r->residual=norm(errors(*r,next,f))/scale;r->iteration++;
        bool converged=std::abs(r->de)<r->etol&&r->drms<r->dtol&&r->residual<r->rtol;
        if(converged){
            Pair canonical=diagonalize(*r,f),canonical_f=fock(*r,canonical);
            const double canonical_e=energy(*r,canonical,canonical_f),canonical_residual=norm(errors(*r,canonical,canonical_f))/scale;
            converged=difference(canonical,next)/scale<r->dtol&&difference(canonical,r->d)/scale<r->dtol&&std::abs(canonical_e-r->energy)<r->etol&&std::abs(canonical_e-e)<r->etol&&canonical_residual<r->rtol;
            next=std::move(canonical);f=std::move(canonical_f);e=canonical_e;r->de=e-r->energy;r->drms=difference(next,r->d)/scale;r->residual=canonical_residual;
        }
        if(!std::isfinite(e)||!next[0].allFinite()||!next[1].allFinite()||!f[0].allFinite()||!f[1].allFinite())throw std::runtime_error("UHF produced a nonfinite density or energy.");
        r->d=std::move(next);r->f=std::move(f);r->energy=e;r->status=converged?1:r->iteration>=r->maxiter?2:0;return r->status;
    }catch(const std::exception &e){last_error=e.what();r->status=-1;return -1;}
}
extern "C" int qc_uhf_nao(const QCUHF *r){return r?r->n:0;}
extern "C" int qc_uhf_status(const QCUHF *r){return r?r->status:-1;}
extern "C" int qc_uhf_iteration(const QCUHF *r){return r?r->iteration:0;}
extern "C" double qc_uhf_scalar(const QCUHF *r,int kind){
    if(!r)return std::numeric_limits<double>::quiet_NaN();
    switch(kind){case 0:return r->energy;case 1:return r->de;case 2:return r->drms;case 3:return r->residual;case 4:return r->overlap_min;
        case 5:{const Matrix overlap=r->c[0].leftCols(r->occ[0]).transpose()*r->s*r->c[1].leftCols(r->occ[1]);const double sz=.5*(r->occ[0]-r->occ[1]);return std::max(0.0,sz*(sz+1)+r->occ[1]-overlap.squaredNorm());}
        case 6:{const double s=.5*(r->occ[0]-r->occ[1]);return s*(s+1);}case 7:return r->mixed?1:0;default:return std::numeric_limits<double>::quiet_NaN();}
}
extern "C" const double *qc_uhf_data(const QCUHF *r,int kind){
    if(!r)return nullptr;
    switch(kind){case 0:return r->d[0].data();case 1:return r->d[1].data();case 2:return r->f[0].data();case 3:return r->f[1].data();case 4:return r->c[0].data();case 5:return r->c[1].data();case 6:return r->eps[0].data();case 7:return r->eps[1].data();case 8:return r->occupations[0].data();case 9:return r->occupations[1].data();default:return nullptr;}
}
extern "C" void qc_uhf_free(QCUHF *r){delete r;}

extern "C" QCUHF *qc_uhf_create(const QCResult *i,int na,int nb,int mi,int di,int g,double e,double d,double r){return qc_uhf_create_df(i,nullptr,na,nb,mi,di,g,e,d,r);}

extern "C" int qc_uhf_df_rank(const QCUHF *r){return r&&r->df?r->df->rank:0;}
