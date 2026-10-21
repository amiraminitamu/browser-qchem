// SPDX-License-Identifier: AGPL-3.0-only
#ifndef LORENTZIAN_DF_FACTORS_H
#define LORENTZIAN_DF_FACTORS_H
#include "df_integrals.h"
#include <Eigen/Dense>
#include <stdexcept>
#include <algorithm>
using DFMatrix=Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>;
struct DFFactors {
  DFMatrix b; int n,rank; double minimum;
  explicit DFFactors(const QCDF *df,int expected) {
    if(!df||qc_df_nao(df)!=expected)throw std::runtime_error("DF and reference orbital dimensions differ.");
    n=expected;const int m=qc_df_naux(df);
    const Eigen::Map<const DFMatrix> metric(qc_df_data(df,0),m,m);
    Eigen::SelfAdjointEigenSolver<DFMatrix> eig((.5*(metric+metric.transpose())).eval());
    if(eig.info()!=Eigen::Success||!eig.eigenvalues().allFinite())throw std::runtime_error("DF Coulomb metric diagonalization failed.");
    minimum=eig.eigenvalues()[0];
    if(minimum < -std::max(1e-8,1e-10*eig.eigenvalues()[m-1]))throw std::runtime_error("DF Coulomb metric is not positive semidefinite.");
    rank=(eig.eigenvalues().array()>1e-10).count();
    if(!rank)throw std::runtime_error("All auxiliary metric directions were removed.");
    DFMatrix w=eig.eigenvectors().rightCols(rank);
    for(int p=0;p<rank;++p)w.col(p)/=std::sqrt(eig.eigenvalues()[m-rank+p]);
    b=Eigen::Map<const DFMatrix>(qc_df_data(df,1),n*n,m)*w;
    if(!b.allFinite())throw std::runtime_error("DF whitening produced nonfinite factors.");
  }
  DFMatrix coulomb(const DFMatrix &d) const {
    Eigen::VectorXd q=b.transpose()*Eigen::Map<const Eigen::VectorXd>(d.data(),n*n);
    Eigen::VectorXd j=b*q;
    return Eigen::Map<const DFMatrix>(j.data(),n,n);
  }
  DFMatrix exchange(const DFMatrix &occupied) const {
    // L(mu,i,Q)=sum_nu B(mu,nu,Q) C(nu,i), K=L L^T.
    DFMatrix l(n,occupied.cols()*rank);
    for(int mu=0;mu<n;++mu){DFMatrix block=occupied.transpose()*b.middleRows(mu*n,n);l.row(mu)=Eigen::Map<const Eigen::RowVectorXd>(block.data(),block.size());}
    return l*l.transpose();
  }
};
#endif
