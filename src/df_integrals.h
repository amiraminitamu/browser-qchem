/* SPDX-License-Identifier: AGPL-3.0-only */
#ifndef LORENTZIAN_DF_INTEGRALS_H
#define LORENTZIAN_DF_INTEGRALS_H
#ifdef __cplusplus
extern "C" {
#endif
typedef struct QCDF QCDF;
/* Spherical real GTOs. Orbital shells l<=2; auxiliary shells l<=6.
 * Shell records, coordinates, coefficients have the same convention as qc_compute.
 * Evaluates genuine (mu nu|P) and (P|Q), independently of four-center ERIs. */
QCDF *qc_df_compute(int natm,const int *charges,const double *xyz,
    int nshell,const int *shells,int nprim,const double *exponents,const double *coefficients,
    int nauxshell,const int *auxshells,int nauxprim,const double *auxexponents,const double *auxcoefficients);
const char *qc_df_error(void);
int qc_df_nao(const QCDF *);
int qc_df_naux(const QCDF *);
/* kind 0: metric J[P,Q] (M*M); kind 1: (mu nu|P) [(mu*N+nu)*M+P]. */
const double *qc_df_data(const QCDF *,int kind);
void qc_df_free(QCDF *);
#ifdef __cplusplus
}
#endif
#endif
