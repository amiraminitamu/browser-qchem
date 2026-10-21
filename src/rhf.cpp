// SPDX-License-Identifier: AGPL-3.0-only
#include "rhf.h"
#include "df_factors.h"
#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

using Matrix = Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>;
using Vector = Eigen::VectorXd;
static std::string last_error;

struct QCRHF {
    int n, nocc, maxiter, iteration=0, status=0;
    bool use_diis;
    double enuc, energy, de=0, drms=0, residual=0, overlap_min;
    double etol, dtol, rtol;
    std::unique_ptr<DFFactors> df;
    const double *eri;
    const QCResult *integrals;
    Matrix s, h, x, d, f, c;
    Vector eps, occupations;
    std::vector<Matrix> fock_history, error_history;
};

static Matrix build_fock(const QCRHF &r, const Matrix &d)
{
    Matrix f=r.h;
    if(r.df) return (r.h+r.df->coulomb(d)-r.df->exchange(r.c.leftCols(r.nocc))).eval();
    const size_t n=r.n;
    for (size_t i=0; i<n; ++i) for (size_t j=0; j<n; ++j) {
        double g=0;
        for (size_t k=0; k<n; ++k) for (size_t l=0; l<n; ++l)
            g += d(k,l)*(r.eri[((i*n+j)*n+k)*n+l] - 0.5*r.eri[((i*n+k)*n+j)*n+l]);
        f(i,j)+=g;
    }
    return (0.5*(f+f.transpose())).eval();
}
static double total_energy(const QCRHF &r, const Matrix &d, const Matrix &f)
{
    return 0.5*(d.array()*(r.h+f).array()).sum()+r.enuc;
}
static Matrix commutator(const QCRHF &r, const Matrix &d, const Matrix &f)
{
    return r.x.transpose()*(f*d*r.s-r.s*d*f)*r.x;
}
static Matrix diagonalize(QCRHF &r, const Matrix &f)
{
    Matrix transformed=r.x.transpose()*f*r.x;
    transformed=(0.5*(transformed+transformed.transpose())).eval();
    Eigen::SelfAdjointEigenSolver<Matrix> eig(transformed);
    if (eig.info()!=Eigen::Success) throw std::runtime_error("Fock diagonalization failed.");
    r.c=r.x*eig.eigenvectors(); r.eps=eig.eigenvalues();
    Matrix occupied=r.c.leftCols(r.nocc);
    return 2.0*occupied*occupied.transpose();
}
static Matrix diis_fock(QCRHF &r)
{
    if (!r.use_diis) return r.f;
    r.fock_history.push_back(r.f);
    r.error_history.push_back(commutator(r,r.d,r.f));
    if (r.fock_history.size()>8) {
        r.fock_history.erase(r.fock_history.begin());
        r.error_history.erase(r.error_history.begin());
    }
    int m=static_cast<int>(r.fock_history.size());
    if (m<2) return r.f;
    Matrix b=Matrix::Zero(m+1,m+1);
    for (int i=0;i<m;++i) for (int j=0;j<m;++j)
        b(i,j)=(r.error_history[i].array()*r.error_history[j].array()).sum();
    double scale=b.topLeftCorner(m,m).cwiseAbs().maxCoeff();
    if (scale<1e-28) return r.f;
    b.topLeftCorner(m,m)/=scale;
    b.block(0,m,m,1).setConstant(-1); b.block(m,0,1,m).setConstant(-1);
    Vector rhs=Vector::Zero(m+1); rhs(m)=-1;
    Eigen::FullPivLU<Matrix> solve(b);
    solve.setThreshold(1e-12);
    if (!solve.isInvertible()) return r.f;
    Vector coefficients=solve.solve(rhs);
    if (!coefficients.allFinite() || coefficients.head(m).cwiseAbs().maxCoeff()>1e6 ||
        (b*coefficients-rhs).norm()>1e-7) return r.f;
    Matrix f=Matrix::Zero(r.n,r.n);
    for (int i=0;i<m;++i) f+=coefficients(i)*r.fock_history[i];
    return f;
}

extern "C" const char *qc_rhf_error(void) { return last_error.c_str(); }
extern "C" QCRHF *qc_rhf_create_df(const QCResult *integrals,const QCDF *df,int nelectron,int maxiter,
                                int use_diis,double etol,double dtol,double rtol)
{
    last_error.clear();
    try {
        const int n=qc_nao(integrals);
        if (n<1 || n>40 || nelectron<2 || nelectron%2 || nelectron>2*n)
            throw std::runtime_error("RHF requires a positive even electron count that fits the basis (closed-shell singlet).");
        if (maxiter<1 || maxiter>500 || !std::isfinite(etol) || !std::isfinite(dtol) ||
            !std::isfinite(rtol) || etol<=0 || dtol<=0 || rtol<=0)
            throw std::runtime_error("Invalid SCF iteration limit or convergence threshold.");
        auto r=std::make_unique<QCRHF>();
        r->n=n; r->nocc=nelectron/2; r->maxiter=maxiter; r->use_diis=use_diis!=0;
        r->etol=etol; r->dtol=dtol; r->rtol=rtol; r->enuc=qc_enuc(integrals);
        r->eri=qc_data(integrals,3); r->integrals=integrals;
        if(df)r->df=std::make_unique<DFFactors>(df,n);
        else if(!r->eri)throw std::runtime_error("Exact SCF requires four-center integrals.");
        r->s=Eigen::Map<const Matrix>(qc_data(integrals,0),n,n);
        r->h=Eigen::Map<const Matrix>(qc_data(integrals,1),n,n)
             +Eigen::Map<const Matrix>(qc_data(integrals,2),n,n);
        Eigen::SelfAdjointEigenSolver<Matrix> eig(r->s);
        if (eig.info()!=Eigen::Success) throw std::runtime_error("Overlap diagonalization failed.");
        r->overlap_min=eig.eigenvalues().minCoeff();
        if (r->overlap_min<1e-8)
            throw std::runtime_error("Overlap is nearly linearly dependent (smallest eigenvalue < 1e-8). Adjust the geometry.");
        r->x=eig.eigenvectors()*eig.eigenvalues().cwiseSqrt().cwiseInverse().asDiagonal()*eig.eigenvectors().transpose();
        r->d=diagonalize(*r,r->h);
        r->f=build_fock(*r,r->d);
        r->energy=total_energy(*r,r->d,r->f);
        r->residual=commutator(*r,r->d,r->f).norm()/n;
        r->occupations=Vector::Zero(n); r->occupations.head(r->nocc).setConstant(2);
        return r.release();
    } catch (const std::exception &error) { last_error=error.what(); return nullptr; }
}
extern "C" int qc_rhf_step(QCRHF *r)
{
    if (!r) { last_error="Missing RHF handle."; return -1; }
    if (r->status!=0) return r->status;
    try {
        Matrix next=diagonalize(*r,diis_fock(*r));
        Matrix f=build_fock(*r,next);
        double energy=total_energy(*r,next,f);
        r->de=energy-r->energy;
        r->drms=(next-r->d).norm()/r->n;
        r->residual=commutator(*r,next,f).norm()/r->n;
        r->iteration++;
        bool converged=std::abs(r->de)<r->etol && r->drms<r->dtol && r->residual<r->rtol;
        if (converged) {
            // Canonicalize using the physical Fock and verify a final pure SCF
            // step. Never report orbitals from an unchecked extrapolated Fock.
            Matrix canonical=diagonalize(*r,f);
            Matrix canonical_f=build_fock(*r,canonical);
            double canonical_energy=total_energy(*r,canonical,canonical_f);
            double canonical_residual=commutator(*r,canonical,canonical_f).norm()/r->n;
            converged=(canonical-next).norm()/r->n<r->dtol &&
                (canonical-r->d).norm()/r->n<r->dtol &&
                std::abs(canonical_energy-r->energy)<r->etol &&
                std::abs(canonical_energy-energy)<r->etol && canonical_residual<r->rtol;
            next=std::move(canonical); f=std::move(canonical_f); energy=canonical_energy;
            r->de=energy-r->energy; r->drms=(next-r->d).norm()/r->n;
            r->residual=canonical_residual;
        }
        if (!std::isfinite(energy) || !next.allFinite() || !f.allFinite())
            throw std::runtime_error("SCF produced a nonfinite density or energy.");
        r->d=std::move(next); r->f=std::move(f); r->energy=energy;
        r->status=converged ? 1 : (r->iteration>=r->maxiter ? 2 : 0);
        return r->status;
    } catch (const std::exception &error) { last_error=error.what(); r->status=-1; return -1; }
}
extern "C" int qc_rhf_status(const QCRHF *r) { return r ? r->status : -1; }
extern "C" int qc_rhf_iteration(const QCRHF *r) { return r ? r->iteration : 0; }
extern "C" double qc_rhf_scalar(const QCRHF *r,int kind)
{
    if (!r) return std::numeric_limits<double>::quiet_NaN();
    switch (kind) { case 0:return r->energy; case 1:return r->de; case 2:return r->drms;
        case 3:return r->residual; case 4:return r->overlap_min; default:return std::numeric_limits<double>::quiet_NaN(); }
}
extern "C" const double *qc_rhf_data(const QCRHF *r,int kind)
{
    if (!r) return nullptr;
    switch (kind) { case 0:return r->d.data(); case 1:return r->f.data(); case 2:return r->c.data();
        case 3:return r->eps.data(); case 4:return r->occupations.data(); default:return nullptr; }
}
extern "C" void qc_rhf_free(QCRHF *r) { delete r; }


extern "C" const QCResult *qc_rhf_integrals(const QCRHF *r) { return r ? r->integrals : nullptr; }

extern "C" QCRHF *qc_rhf_create(const QCResult *i,int ne,int mi,int di,double e,double d,double r){return qc_rhf_create_df(i,nullptr,ne,mi,di,e,d,r);}

extern "C" int qc_rhf_df_rank(const QCRHF *r){return r&&r->df?r->df->rank:0;}
