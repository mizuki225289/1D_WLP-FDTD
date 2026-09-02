#include "const_region.h"
#include "head.h"
#include "PMLparameter.h"
#include "r_theta.h"

constexpr double b0 {2.0 / MU0 / s};

void calc_Hphi(double* Hphi_q, double* Er_q, double* Hphi_sum, double* hphi_sum) {
    /*非PML*/
    for(int j=0; j < N_THETA - L; j++) {
        Hphi_q[j] = b0 / EARTH_R / Delta_THETA * (Er_q[j+1] - Er_q[j]) - 2.0 * Hphi_sum[j];
    }

    /*PML*/
    for(int j = N_THETA - L; j < N_THETA; j++) {
        Hphi_q[j]
            = xi_theta(j+0.5) * b0 / EARTH_R / Delta_THETA * (Er_q[j+1] - Er_q[j])
            + rho_theta(j+0.5) * hphi_sum[j]
            - tau_theta(j+0.5) * Hphi_sum[j];
    }
}