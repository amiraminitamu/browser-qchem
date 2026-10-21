/* SPDX-License-Identifier: AGPL-3.0-only */
#ifndef BROWSER_QCHEM_UHF_H
#define BROWSER_QCHEM_UHF_H
#include "integrals.h"
#include "df_integrals.h"
#ifdef __cplusplus
extern "C" {
#endif
typedef struct QCUHF QCUHF;
/* Borrow integrals until this handle is freed. nalpha >= nbeta >= 0.
 * guess: 0 core, 1 opposite HOMO/LUMO rotations for equal spin populations.
 * Status: 0 iterating, 1 converged, 2 iteration limit, -1 error. */
QCUHF *qc_uhf_create(const QCResult *,int nalpha,int nbeta,int maxiter,int diis,
                     int guess,double etol,double dtol,double rtol);
/* DF factors are copied at construction; the QCDF may then be freed. */
QCUHF *qc_uhf_create_df(const QCResult *,const QCDF *,int na,int nb,int maxiter,int diis,int guess,double etol,double dtol,double rtol);
int qc_uhf_df_rank(const QCUHF *);
const char *qc_uhf_error(void);
int qc_uhf_step(QCUHF *);
int qc_uhf_nao(const QCUHF *);
int qc_uhf_status(const QCUHF *);
int qc_uhf_iteration(const QCUHF *);
/* 0 E, 1 delta E, 2 spin-averaged density RMS, 3 spin-averaged residual RMS,
 * 4 min overlap eigenvalue, 5 <S^2>, 6 target S(S+1), 7 mixed guess applied. */
double qc_uhf_scalar(const QCUHF *,int kind);
/* Row-major: 0 Da,1 Db,2 Fa,3 Fb,4 Ca,5 Cb,6 epsa,7 epsb,8 occa,9 occb.
 * Matrices have n*n entries; orbital vectors have n. Spin occupations are 1/0. */
const double *qc_uhf_data(const QCUHF *,int kind);
void qc_uhf_free(QCUHF *);
#ifdef __cplusplus
}
#endif
#endif
