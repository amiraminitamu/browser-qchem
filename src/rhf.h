/* SPDX-License-Identifier: AGPL-3.0-only */
#ifndef BROWSER_QCHEM_RHF_H
#define BROWSER_QCHEM_RHF_H
#include "integrals.h"
#include "df_integrals.h"
#ifdef __cplusplus
extern "C" {
#endif

typedef struct QCRHF QCRHF;
/* The integral handle must outlive the SCF handle. D is spin-summed, with
 * occupation 2 per occupied spatial orbital. Status: 0 iterating, 1 converged,
 * 2 max iterations, -1 error.
 * One call to qc_rhf_step performs one SCF iteration; inspect progress in between.
 */
QCRHF *qc_rhf_create(const QCResult *integrals, int nelectron, int maxiter,
                     int use_diis, double energy_tol, double density_tol,
                     double residual_tol);
/* DF factors are copied at construction; the QCDF may then be freed. */
QCRHF *qc_rhf_create_df(const QCResult *,const QCDF *,int nelectron,int maxiter,int diis,double etol,double dtol,double rtol);
int qc_rhf_df_rank(const QCRHF *);
const char *qc_rhf_error(void);
int qc_rhf_step(QCRHF *rhf);
int qc_rhf_status(const QCRHF *rhf);
int qc_rhf_iteration(const QCRHF *rhf);
/* Scalars: 0 total energy, 1 energy change, 2 density RMS change,
 * 3 orthonormal-basis commutator RMS, 4 smallest overlap eigenvalue.
 */
double qc_rhf_scalar(const QCRHF *rhf, int kind);
/* Arrays, row-major: 0 D, 1 physical F(D), 2 MO coefficients C[AO,MO],
 * 3 orbital energies (n), 4 occupations (n).
 */
const double *qc_rhf_data(const QCRHF *rhf, int kind);
const QCResult *qc_rhf_integrals(const QCRHF *rhf);
void qc_rhf_free(QCRHF *rhf);
#ifdef __cplusplus
}
#endif
#endif

