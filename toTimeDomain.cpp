#include "const_sq.h"
#include "const_time.h"

void toTimeDomain(double* timeDomain, double* WLPcoef, double** Laguerre, int Nt) {
    for(int t=0; t < Nt; t++) {
        double tmp = 0;
        for(int p=0; p < q_MAX; p++) {
            tmp += WLPcoef[p] * Laguerre[p][t];
        }
        timeDomain[t] = tmp;
    }
}

void toTimeDomain_Er(int q, double** Er_t, double* Er_q, double** Laguerre) {
    for(int j=0; j <= N_THETA; j++) {
        for(int t=0; t < Nt_output; t++) {
            Er_t[j][t] += Er_q[j] * Laguerre[q][t];
        }
    }
}