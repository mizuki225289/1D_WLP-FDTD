#include "const_region.h"

void update_sum_r(double* Er_q, double* Er_sum) {
    for(int j=0 ; j < N_THETA + 1; j++) {
        Er_sum[j] += Er_q[j];
    }
}

void update_sum_phi(double* Hphi_q, double* Hphi_sum) {
    for(int j=0; j < N_THETA; j++) {
        Hphi_sum[j] += Hphi_q[j];
    }
}