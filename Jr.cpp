#include "const.h"
#include "const_region.h"
#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>

/*印加電流について*/
constexpr double f = 22.2e3;
constexpr double omega0 = 2.0 * M_PI * f;
constexpr double sigma = 1.0 / f;
constexpr double tau0 = 6.0 * sigma;

double Jr(double t) {
    return (t - tau0) * exp(-(t - tau0) * (t - tau0) / (2.0 * sigma * sigma)) / sigma;
}