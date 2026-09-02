#ifndef __PML_PARAMETER__
#define __PML_PARAMETER__

#define _USE_MATH_DEFINES
#include <math.h>
#include <iostream>
#include "head.h"
#include "const.h"
#include "const_region.h"
#include "const_sq.h"
#include "r_theta.h"

constexpr double M {1.0};
constexpr double R {1.0e-6};

const double sigma_theta_max = -(M + 1.0) * C0 * log(R) / 2.0 / L / EARTH_R / Delta_THETA;

template <typename ValueType>
double sigma_theta(ValueType j) {
    if(j < N_THETA - L) {
        return 0;
    } else {
        double tmp = (1.0 * j - N_THETA + L) / L;
        return sigma_theta_max * pow(tmp, M); 
    }
}

template <typename ValueType>
double beta_theta(ValueType j) {
    return 1.0 / (sigma_theta(j) + s / 2.0);
}

template <typename ValueType>
double gamma_theta(ValueType j) {
    return (sigma_theta(j) + s / 2.0);
}

template <typename ValueType>
double kappa_theta(ValueType j) {
    return 1;
}

template <typename ValueType>
double alpha_theta(ValueType j) {
    return 0;
}

template <typename ValueType>
double zeta_theta(ValueType j) {
    return alpha_theta(j) * kappa_theta(j) + sigma_theta(j) + kappa_theta(j) * s / 2.0;
}

template <typename ValueType>
double eta_theta(ValueType j) {
    return s / 2.0 + alpha_theta(j);
}

template <typename ValueType>
double xi_theta(ValueType j) {
    return eta_theta(j) / zeta_theta(j);
}

template <typename ValueType>
double rho_theta(ValueType j) {
    return (s - 2.0 * eta_theta(j)) / zeta_theta(j);
}

template <typename ValueType>
double tau_theta(ValueType j) {
    return kappa_theta(j) * s / zeta_theta(j);
}

#endif