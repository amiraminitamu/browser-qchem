/* SPDX-License-Identifier: AGPL-3.0-only */
#ifndef BROWSER_QCHEM_INTEGRALS_H
#define BROWSER_QCHEM_INTEGRALS_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct QCResult QCResult;

/* Coordinates in bohr. Each shell record contains [atom_index, l, nprim].
 * One contraction per shell, l=0..2, normalized primitives as input convention.
 * Concatenate exponents and contraction coefficients in shell order.
 * The caller owns input buffers; the returned handle owns all result buffers.
 */
QCResult *qc_compute(int natm, const int *charges, const double *xyz,
                     int nshell, const int *shells, int nprim,
                     const double *exponents, const double *coefficients);
/* dense=0 computes only S/T/V; preserves the libcint context for on-demand ERIs. */
QCResult *qc_compute_mode(int dense,int natm,const int *charges,const double *xyz,int nshell,const int *shells,int nprim,const double *exponents,const double *coefficients);
double qc_eri(const QCResult *,int i,int j,int k,int l);
int qc_eri_transform(const QCResult *,const double *coefficient_column,double *out);
const char *qc_error(void);
int qc_nao(const QCResult *result);
double qc_enuc(const QCResult *result);
/* kind=0 overlap, 1 kinetic, 2 nuclear attraction, 3 ERI.
 * All arrays are row-major. ERI is chemists' (ij|kl), kets are real spherical.
 */
const double *qc_data(const QCResult *result, int kind);
void qc_free(QCResult *result);

#ifdef __cplusplus
}
#endif

#endif
