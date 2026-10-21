// SPDX-License-Identifier: AGPL-3.0-only
#include "mp2.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>
struct QCMP2 {
    int n,o,v,iteration=0,status=0;
    const double *eri,*c,*eps;
    const QCResult *integrals;
    double reference,os=0,ss=0,min_denominator=std::numeric_limits<double>::infinity(),gap;
};
static std::string last_error;
extern "C" QCMP2 *qc_mp2_create(const QCRHF *reference) {
    last_error.clear();
    try {
        if (!reference || qc_rhf_status(reference)!=1) throw std::runtime_error("MP2 requires a converged RHF reference.");
        auto m=std::make_unique<QCMP2>();
        const QCResult *integrals=qc_rhf_integrals(reference);
        m->n=qc_nao(integrals); m->o=0;
        const double *occupation=qc_rhf_data(reference,4);
        for (int i=0;i<m->n;++i) if (occupation[i]>1) ++m->o;
        m->integrals=integrals; m->v=m->n-m->o; m->eri=qc_data(integrals,3);
        m->c=qc_rhf_data(reference,2); m->eps=qc_rhf_data(reference,3);
        m->reference=qc_rhf_scalar(reference,0);
        m->gap=m->v ? m->eps[m->o]-m->eps[m->o-1] : std::numeric_limits<double>::quiet_NaN();
        if (!m->v) m->status=1;
        return m.release();
    } catch (const std::exception &e) { last_error=e.what(); return nullptr; }
}
extern "C" int qc_mp2_step(QCMP2 *m) {
    if (!m) { last_error="Invalid MP2 handle."; return -1; }
    if (m->status) return m->status;
    try {
        const int n=m->n,o=m->o,v=m->v,i=m->iteration;
        const double *c=m->c;
        // Four quarter-transformations: (mu nu|lambda sigma) ->
        // (i nu|lambda sigma) -> (i a|lambda sigma) -> (i a|j sigma) -> (i a|j b).
        // Keep i fixed, so the largest extra array is only n^3.
        std::vector<double> x(n*n*n,0),y(v*n*n,0),z(v*o*n,0),g(v*o*v,0);
        if(m->eri) {
        for(int mu=0;mu<n;++mu) for(int nu=0;nu<n;++nu)
            for(int la=0;la<n;++la) for(int si=0;si<n;++si)
                x[(nu*n+la)*n+si]+=c[mu*n+i]*m->eri[((mu*n+nu)*n+la)*n+si];
        } else {
            std::vector<double> column(n);for(int mu=0;mu<n;++mu)column[mu]=c[mu*n+i];
            if(!qc_eri_transform(m->integrals,column.data(),x.data()))throw std::runtime_error("Direct MP2 integral transformation failed.");
        }
        for(int a=0;a<v;++a) for(int nu=0;nu<n;++nu)
            for(int la=0;la<n;++la) for(int si=0;si<n;++si)
                y[(a*n+la)*n+si]+=c[nu*n+o+a]*x[(nu*n+la)*n+si];
        for(int a=0;a<v;++a) for(int j=0;j<o;++j)
            for(int la=0;la<n;++la) for(int si=0;si<n;++si)
                z[(a*o+j)*n+si]+=c[la*n+j]*y[(a*n+la)*n+si];
        for(int a=0;a<v;++a) for(int j=0;j<o;++j)
            for(int b=0;b<v;++b) for(int si=0;si<n;++si)
                g[(a*o+j)*v+b]+=c[si*n+o+b]*z[(a*o+j)*n+si];
        double os=0,ss=0,min_denominator=m->min_denominator;
        for(int j=0;j<o;++j) for(int a=0;a<v;++a) for(int b=0;b<v;++b) {
            const double d=m->eps[i]+m->eps[j]-m->eps[o+a]-m->eps[o+b];
            if (!std::isfinite(d) || d>=-1e-10) throw std::runtime_error("MP2 denominator is non-negative or too small. The RHF reference is unsuitable for MP2.");
            min_denominator=std::min(min_denominator,-d);
            const double direct=g[(a*o+j)*v+b],exchange=g[(b*o+j)*v+a];
            os+=direct*direct/d; ss+=direct*(direct-exchange)/d;
        }
        if(!std::isfinite(os+ss)) throw std::runtime_error("MP2 correlation energy is not finite.");
        m->os+=os;m->ss+=ss;m->min_denominator=min_denominator;
        if(++m->iteration==o) m->status=1;
    } catch (const std::exception &e) {last_error=e.what();m->status=-1;}
    return m->status;
}
extern "C" const char *qc_mp2_error(void) {return last_error.c_str();}
extern "C" int qc_mp2_status(const QCMP2 *m) {return m?m->status:-1;}
extern "C" int qc_mp2_iteration(const QCMP2 *m) {return m?m->iteration:0;}
extern "C" int qc_mp2_nocc(const QCMP2 *m) {return m?m->o:0;}
extern "C" int qc_mp2_nvirt(const QCMP2 *m) {return m?m->v:0;}
extern "C" double qc_mp2_scalar(const QCMP2 *m,int kind) {
    if(m) switch(kind) {case 0:return m->reference+m->os+m->ss;case 1:return m->os+m->ss;
        case 2:return m->os;case 3:return m->ss;case 4:return m->min_denominator;case 5:return m->gap;case 6:return m->reference;}
    return std::numeric_limits<double>::quiet_NaN();
}
extern "C" void qc_mp2_free(QCMP2 *m) {delete m;}
