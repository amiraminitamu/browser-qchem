/* SPDX-License-Identifier: AGPL-3.0-only */
#ifndef LORENTZIAN_DF_MP2_H
#define LORENTZIAN_DF_MP2_H
#include "df_integrals.h"
#include "rhf.h"
#include "uhf.h"
#ifdef __cplusplus
extern "C" {
#endif
typedef struct QCDFMP2 QCDFMP2;
/* Canonical, all-electron Coulomb-metric RI-MP2. Reference must be converged.
 * Only the construction borrows reference/DF data; the returned handle owns its
 * transformed three-index factors and orbital energies. No amplitudes or
 * four-index MO tensor are stored. Absolute metric eigenvalue cutoff: 1e-10.
 * Status: 0 running, 1 complete, -1 error. Each step contracts an occupied block. */
QCDFMP2 *qc_dfmp2_rhf_create(const QCRHF *,const QCDF *);
QCDFMP2 *qc_dfmp2_uhf_create(const QCUHF *,const QCDF *);
const char *qc_dfmp2_error(void);
int qc_dfmp2_step(QCDFMP2 *);
int qc_dfmp2_status(const QCDFMP2 *);
int qc_dfmp2_iteration(const QCDFMP2 *);
int qc_dfmp2_steps(const QCDFMP2 *);
int qc_dfmp2_nocc(const QCDFMP2 *,int spin);
int qc_dfmp2_nvirt(const QCDFMP2 *,int spin);
int qc_dfmp2_rank(const QCDFMP2 *);
/* 0 total,1 correlation,2 opposite-spin,3 same-spin,4 minimum |denominator|,
 * 5 minimum reference spin gap,6 reference energy,7 alpha-alpha,8 beta-beta,
 * 9 metric cutoff,10 smallest metric eigenvalue. */
double qc_dfmp2_scalar(const QCDFMP2 *,int kind);
void qc_dfmp2_free(QCDFMP2 *);
#ifdef __cplusplus
}
#endif
#endif
