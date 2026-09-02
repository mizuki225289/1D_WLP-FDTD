#include "const.h"
#include "const_sq.h"
#include "const_region.h"
#include "PMLParameter.h"
#include "head.h"
#include "r_theta.h"

/*非磁化プラズマ*/
constexpr double a0 {2.0 / EPS0 / s};
constexpr double b0 {2.0 / MU0 / s};
constexpr double xi0 {a0 * b0};

void update_beta(double* b, double* Hphi_sum, double* Er_sum, double* hphi_sum, double* er_sum) {
    /*左端　j=0 回転中心　PMLなし*/
    b[0] = -16.0 / EPS0 / s / EARTH_R / Delta_THETA * Hphi_sum[0] - 2.0 * Er_sum[0];

    /*通常領域*/
    for(int j=1; j < N_THETA - L; j++) {
        b[j]
        = -xi0 / EARTH_R / EARTH_R / sin_j(j) / Delta_THETA * (
            Hphi_sum[j] - Hphi_sum[j-1]
        )
        - 2.0 * Er_sum[j];
    }

    /*PML領域*/
    for(int j = N_THETA - L; j < N_THETA; j++) {
        b[j] 
        = xi_theta(j) * a0 / EARTH_R / sin_j(j) / Delta_THETA * (
            + sin_j(j+0.5) * rho_theta(j+0.5) * hphi_sum[j]
            - sin_j(j-0.5) * rho_theta(j-0.5) * hphi_sum[j-1]
        ) - xi_theta(j) * a0 / EARTH_R / sin_j(j) / Delta_THETA * (
            + sin_j(j+0.5) * tau_theta(j+0.5) * Hphi_sum[j]
            - sin_j(j-0.5) * tau_theta(j-0.5) * Hphi_sum[j-1]
        )
        + rho_theta(j) * er_sum[j]
        - tau_theta(j) * Er_sum[j];
    }

    /*右端 PEC j = N_THETA*/
    b[N_THETA] = 0;
}