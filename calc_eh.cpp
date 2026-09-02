#include "head.h"
#include "const_region.h"
#include "PMLparameter.h"
#include "r_theta.h"

constexpr double a0 {2.0 / EPS0 / s};
constexpr double b0 {2.0 / MU0 / s};

void calc_er_q(double* er_q, double* er_sum, double* Er_q, double* Er_sum) {
    for(int j=0; j < N_THETA + 1; j++) {
        er_q[j] 
            = (zeta_theta(j) / eta_theta(j)) * Er_q[j]
            + (kappa_theta(j) * s / eta_theta(j)) * Er_sum[j]
            - (s / eta_theta(j)) * er_sum[j];
    }
}

void calc_hphi_q(double* hphi_q, double* hphi_sum, double* Hphi_q, double* Hphi_sum) {
    for(int j=0; j < N_THETA; j++) {
        hphi_q[j] 
            = (zeta_theta(j+0.5) / eta_theta(j+0.5)) * Hphi_q[j]
            + (kappa_theta(j+0.5) * s / eta_theta(j+0.5)) * Hphi_sum[j]
            - (s / eta_theta(j+0.5)) * hphi_sum[j];
    }
}