/* SPDX-License-Identifier: AGPL-3.0-only */
#include <math.h>
#include <stddef.h>
#include "cint.h"
#include "cint_funcs.h"

/* Overlap between normalized H STO-3G 1s orbitals at separation r_bohr.
 * This wrapper only packs libcint inputs; libcint evaluates the integral.
 * Basis data match PySCF's STO-3G H entry (Hehre, Stewart, Pople, 1969).
 * Returns NaN for a nonfinite or negative separation. Zero is allowed as
 * an orbital-normalization check, not as a physical molecular geometry.
 */
double h2_overlap(double r_bohr)
{
    if (!isfinite(r_bohr) || r_bohr < 0.0) return NAN;

    const double exponents[3] = {3.42525091, 0.62391373, 0.16885540};
    const double coefficients[3] = {0.15432897, 0.53532814, 0.44463454};
    int atm[2 * ATM_SLOTS] = {0};
    int bas[2 * BAS_SLOTS] = {0};
    double env[PTR_ENV_START + 12] = {0.0};
    const int coord = PTR_ENV_START;
    const int exp_ptr = coord + 6;
    const int coeff_ptr = exp_ptr + 3;

    /* Normalize the contraction as well as its individual primitives. */
    double norm2 = 0.0;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            double a = exponents[i], b = exponents[j];
            norm2 += coefficients[i] * coefficients[j]
                   * pow(2.0 * sqrt(a * b) / (a + b), 1.5);
        }
    }
    for (int p = 0; p < 3; ++p) {
        env[exp_ptr + p] = exponents[p];
        /* libcint expects radial normalization in its coefficient array. */
        env[coeff_ptr + p] = coefficients[p] / sqrt(norm2)
                          * CINTgto_norm(0, exponents[p]);
    }
    for (int atom = 0; atom < 2; ++atom) {
        atm[atom * ATM_SLOTS + CHARGE_OF] = 1;
        atm[atom * ATM_SLOTS + PTR_COORD] = coord + 3 * atom;
        bas[atom * BAS_SLOTS + ATOM_OF] = atom;
        bas[atom * BAS_SLOTS + ANG_OF] = 0;
        bas[atom * BAS_SLOTS + NPRIM_OF] = 3;
        bas[atom * BAS_SLOTS + NCTR_OF] = 1;
        bas[atom * BAS_SLOTS + PTR_EXP] = exp_ptr;
        bas[atom * BAS_SLOTS + PTR_COEFF] = coeff_ptr;
    }
    env[coord + 5] = r_bohr;
    int shells[2] = {0, 1};
    double overlap = 0.0;
    int1e_ovlp_sph(&overlap, NULL, shells, atm, 2, bas, 2, env, NULL, NULL);
    return overlap;
}
