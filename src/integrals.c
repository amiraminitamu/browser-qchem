/* SPDX-License-Identifier: AGPL-3.0-only */
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "cint.h"
#include "cint_funcs.h"
#include "integrals.h"
#include "limits.h"

#define MAX_ATOMS 32
#define MAX_SHELLS QC_MAX_SHELLS
#define MAX_AO QC_MAX_AO
#define MAX_PRIMITIVES 4096

struct QCResult {
    int nao;
    double enuc;
    double *data[4];
    int natm, nshell;
    int *atoms, *basis;
    double *environment;
    int offsets[MAX_SHELLS+1];
};
static const char *last_error = "";
const char *qc_error(void) { return last_error; }
int qc_nao(const QCResult *r) { return r ? r->nao : 0; }
double qc_enuc(const QCResult *r) { return r ? r->enuc : NAN; }
const double *qc_data(const QCResult *r, int kind)
{
    return r && kind >= 0 && kind < 4 ? r->data[kind] : NULL;
}
void qc_free(QCResult *r)
{
    if (!r) return;
    for (int i = 0; i < 4; ++i) free(r->data[i]);
    free(r->atoms); free(r->basis); free(r->environment);
    free(r);
}

QCResult *qc_compute_mode(int dense, int natm, const int *charges, const double *xyz,
                     int nshell, const int *shells, int nprim,
                     const double *exponents, const double *coefficients)
{
    last_error = "";
    if ((dense != 0 && dense != 1) || natm < 1 || natm > MAX_ATOMS || nshell < 1 || nshell > MAX_SHELLS ||
        nprim < 1 || nprim > MAX_PRIMITIVES || !charges || !xyz || !shells ||
        !exponents || !coefficients) {
        last_error = "Invalid input sizes or missing arrays.";
        return NULL;
    }
    double enuc = 0.0;
    for (int a = 0; a < natm; ++a) {
        if (charges[a] < 1 || charges[a] > 118) {
            last_error = "Nuclear charges must be integers from 1 to 118.";
            return NULL;
        }
        for (int d = 0; d < 3; ++d) {
            if (!isfinite(xyz[3*a+d]) || fabs(xyz[3*a+d]) > 1e6) {
                last_error = "Coordinates must be finite and within 1e6 bohr.";
                return NULL;
            }
        }
        for (int b = 0; b < a; ++b) {
            double distance2 = 0.0;
            for (int d = 0; d < 3; ++d) {
                double delta = xyz[3*a+d] - xyz[3*b+d];
                distance2 += delta * delta;
            }
            if (distance2 < 1e-16) {
                last_error = "Two nuclei coincide or are closer than 1e-8 bohr.";
                return NULL;
            }
            enuc += charges[a] * charges[b] / sqrt(distance2);
        }
    }
    int offsets[MAX_SHELLS+1] = {0};
    int primitive_count = 0, max_dim = 1;
    for (int s = 0; s < nshell; ++s) {
        int atom = shells[3*s], l = shells[3*s+1], count = shells[3*s+2];
        if (atom < 0 || atom >= natm || l < 0 || l > 2 || count < 1 || count > 64) {
            last_error = "Invalid shell: require an existing atom, l=0..2, and 1..64 primitives.";
            return NULL;
        }
        primitive_count += count;
        offsets[s+1] = offsets[s] + 2*l + 1;
        if (2*l+1 > max_dim) max_dim = 2*l+1;
    }
    if (primitive_count != nprim || offsets[nshell] > MAX_AO) {
        last_error = "Primitive count mismatch or more than 256 AOs.";
        return NULL;
    }
    if(dense && 8.0*pow(offsets[nshell],4)>QC_MAX_DENSE_BYTES){last_error="Exact ERI allocation exceeds 256 MiB. Use DF-JK or a smaller basis.";return NULL;}
    for (int p = 0; p < nprim; ++p) {
        if (!isfinite(exponents[p]) || exponents[p] <= 0 || exponents[p] > 1e8 ||
            !isfinite(coefficients[p]) || fabs(coefficients[p]) > 1e8) {
            last_error = "Invalid primitive exponent or coefficient.";
            return NULL;
        }
    }

    QCResult *r = calloc(1, sizeof(*r));
    int *atm = calloc((size_t)natm * ATM_SLOTS, sizeof(int));
    int *bas = calloc((size_t)nshell * BAS_SLOTS, sizeof(int));
    double *env = calloc(PTR_ENV_START + (size_t)3*natm + (size_t)2*nprim, sizeof(double));
    double *buffer = calloc((size_t)max_dim*max_dim*max_dim*max_dim, sizeof(double));
    CINTOpt *optimizer = NULL;
    if (!r || !atm || !bas || !env || !buffer) goto allocation_failed;
    r->nao = offsets[nshell];
    r->enuc = enuc;
    size_t n = (size_t)r->nao, n2 = n*n;
    for (int kind = 0; kind < (dense ? 4 : 3); ++kind) {
        r->data[kind] = calloc(kind == 3 ? n2*n2 : n2, sizeof(double));
        if (!r->data[kind]) goto allocation_failed;
    }
    for (int a = 0; a < natm; ++a) {
        atm[a*ATM_SLOTS + CHARGE_OF] = charges[a];
        atm[a*ATM_SLOTS + PTR_COORD] = PTR_ENV_START + 3*a;
        atm[a*ATM_SLOTS + NUC_MOD_OF] = POINT_NUC;
        memcpy(env + PTR_ENV_START + 3*a, xyz + 3*a, 3*sizeof(double));
    }
    int primitive = 0, env_offset = PTR_ENV_START + 3*natm;
    for (int s = 0; s < nshell; ++s) {
        int l = shells[3*s+1], count = shells[3*s+2];
        bas[s*BAS_SLOTS + ATOM_OF] = shells[3*s];
        bas[s*BAS_SLOTS + ANG_OF] = l;
        bas[s*BAS_SLOTS + NPRIM_OF] = count;
        bas[s*BAS_SLOTS + NCTR_OF] = 1;
        bas[s*BAS_SLOTS + PTR_EXP] = env_offset;
        bas[s*BAS_SLOTS + PTR_COEFF] = env_offset + count;
        double norm2 = 0.0;
        for (int i = 0; i < count; ++i) {
            for (int j = 0; j < count; ++j) {
                double a = exponents[primitive+i], b = exponents[primitive+j];
                norm2 += coefficients[primitive+i] * coefficients[primitive+j]
                       * pow(2*sqrt(a*b)/(a+b), l+1.5);
            }
        }
        if (!isfinite(norm2) || norm2 < 1e-30) {
            last_error = "A basis contraction has zero or invalid norm.";
            goto failed;
        }
        for (int p = 0; p < count; ++p) {
            double exponent = exponents[primitive+p];
            env[env_offset+p] = exponent;
            env[env_offset+count+p] = coefficients[primitive+p]
                * CINTgto_norm(l, exponent) / sqrt(norm2);
        }
        primitive += count;
        env_offset += 2*count;
    }

    CINTIntegralFunction *one_electron[3] = {int1e_ovlp_sph, int1e_kin_sph, int1e_nuc_sph};
    for (int a = 0; a < nshell; ++a) {
        int da = offsets[a+1]-offsets[a];
        for (int b = 0; b < nshell; ++b) {
            int db = offsets[b+1]-offsets[b];
            int shls[4] = {a,b,0,0};
            for (int kind = 0; kind < 3; ++kind) {
                memset(buffer, 0, (size_t)da*db*sizeof(double));
                one_electron[kind](buffer,NULL,shls,atm,natm,bas,nshell,env,NULL,NULL);
                for (int j = 0; j < db; ++j) {
                    for (int i = 0; i < da; ++i) {
                        r->data[kind][(offsets[a]+i)*n+offsets[b]+j] = buffer[i+da*j];
                    }
                }
            }
        }
    }
    if (dense) {
    int2e_optimizer(&optimizer,atm,natm,bas,nshell,env);
    for (int a = 0; a < nshell; ++a) {
        int da = offsets[a+1]-offsets[a];
        for (int b = 0; b < nshell; ++b) {
            int db = offsets[b+1]-offsets[b];
            for (int c = 0; c < nshell; ++c) {
                int dc = offsets[c+1]-offsets[c];
                for (int d = 0; d < nshell; ++d) {
                    int dd = offsets[d+1]-offsets[d];
                    int shls[4] = {a,b,c,d};
                    memset(buffer, 0, (size_t)da*db*dc*dd*sizeof(double));
                    int2e_sph(buffer,NULL,shls,atm,natm,bas,nshell,env,optimizer,NULL);
                    /* libcint shell blocks have their first AO index fastest.
                     * Export row-major ERI[i,j,k,l] with its last index fastest.
                     */
                    for (int l = 0; l < dd; ++l)
                    for (int k = 0; k < dc; ++k)
                    for (int j = 0; j < db; ++j)
                    for (int i = 0; i < da; ++i) {
                        size_t index = (((offsets[a]+i)*n+offsets[b]+j)*n+offsets[c]+k)*n+offsets[d]+l;
                        r->data[3][index] = buffer[i+da*(j+db*(k+dc*l))];
                    }
                }
            }
        }
    }
    CINTdel_optimizer(&optimizer);
    }
    r->natm=natm; r->nshell=nshell;
    r->atoms=atm; r->basis=bas; r->environment=env;
    memcpy(r->offsets,offsets,sizeof(offsets));
    free(buffer);
    return r;

allocation_failed:
    last_error = "Could not allocate integral buffers.";
failed:
    CINTdel_optimizer(&optimizer);
    free(atm); free(bas); free(env); free(buffer);
    qc_free(r);
    return NULL;
}

QCResult *qc_compute(int natm,const int *charges,const double *xyz,int nshell,const int *shells,int nprim,const double *exponents,const double *coefficients) {
    return qc_compute_mode(1,natm,charges,xyz,nshell,shells,nprim,exponents,coefficients);
}
/* Evaluate exactly one shell quartet, retaining no dense ERI tensor. */
double qc_eri(const QCResult *r,int i,int j,int k,int l) {
    if(!r || i<0 || j<0 || k<0 || l<0 || i>=r->nao || j>=r->nao || k>=r->nao || l>=r->nao) return NAN;
    if(r->data[3])return r->data[3][((i*r->nao+j)*r->nao+k)*r->nao+l];
    int indices[4]={i,j,k,l},shls[4],local[4],dim[4];
    for(int a=0;a<4;++a){int s=0;while(r->offsets[s+1]<=indices[a])++s;shls[a]=s;local[a]=indices[a]-r->offsets[s];dim[a]=r->offsets[s+1]-r->offsets[s];}
    double buffer[625]={0}; /* orbital angular momentum <=2 */
    int2e_sph(buffer,NULL,shls,r->atoms,r->natm,r->basis,r->nshell,r->environment,NULL,NULL);
    return buffer[local[0]+dim[0]*(local[1]+dim[1]*(local[2]+dim[2]*local[3]))];
}
/* First quarter transform for conventional MP2. O(N^3) output, no N^4 storage. */
int qc_eri_transform(const QCResult *r,const double *c,double *out) {
    if(!r||!c||!out)return 0;
    const int n=r->nao;memset(out,0,(size_t)n*n*n*sizeof(double));
    double buffer[625];
    for(int a=0;a<r->nshell;++a)for(int b=0;b<r->nshell;++b)
    for(int csh=0;csh<r->nshell;++csh)for(int d=0;d<r->nshell;++d){
        int da=r->offsets[a+1]-r->offsets[a],db=r->offsets[b+1]-r->offsets[b],dc=r->offsets[csh+1]-r->offsets[csh],dd=r->offsets[d+1]-r->offsets[d];
        int shls[4]={a,b,csh,d};memset(buffer,0,(size_t)da*db*dc*dd*sizeof(double));
        int2e_sph(buffer,NULL,shls,r->atoms,r->natm,r->basis,r->nshell,r->environment,NULL,NULL);
        for(int l=0;l<dd;++l)for(int k=0;k<dc;++k)for(int j=0;j<db;++j)for(int i=0;i<da;++i)
            out[((r->offsets[b]+j)*n+r->offsets[csh]+k)*n+r->offsets[d]+l]+=c[r->offsets[a]+i]*buffer[i+da*(j+db*(k+dc*l))];
    }
    return 1;
}
