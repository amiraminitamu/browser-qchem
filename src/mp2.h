/* SPDX-License-Identifier: AGPL-3.0-only */
#ifndef BROWSER_QCHEM_MP2_H
#define BROWSER_QCHEM_MP2_H
#include "rhf.h"
#ifdef __cplusplus
extern "C" {
#endif
typedef struct QCMP2 QCMP2;
/* All-electron conventional MP2. Reference and its integrals must outlive this
 * handle. Each step contracts one occupied orbital. Status 0 running, 1 done,
 * -1 error. No four-index MO tensor or amplitudes are retained. */
QCMP2 *qc_mp2_create(const QCRHF *reference);
const char *qc_mp2_error(void);
int qc_mp2_step(QCMP2 *mp2);
int qc_mp2_status(const QCMP2 *mp2);
int qc_mp2_iteration(const QCMP2 *mp2);
int qc_mp2_nocc(const QCMP2 *mp2);
int qc_mp2_nvirt(const QCMP2 *mp2);
/* Scalars: total, correlation, opposite-spin, same-spin, minimum absolute
 * denominator, reference HOMO-LUMO gap, reference total energy. */
double qc_mp2_scalar(const QCMP2 *mp2,int kind);
void qc_mp2_free(QCMP2 *mp2);
#ifdef __cplusplus
}
#endif
#endif
