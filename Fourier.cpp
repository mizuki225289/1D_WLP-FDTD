#define _USE_MATH_DEFINES
#include <iostream>
#include <fstream>
#include <cmath>
#include <complex>
#include <string>
#include <vector>
#include "head.h"
#include "const_time.h"

std::complex <double> zj {0.0, 1.0};

std::complex <double> Fourier(const double* data, double omega, int N) {
    std::complex <double> result = 0;
    for(int i=0; i < N; i++) {
        result += Delta_t_pulse * data[i] * exp( -zj * double(i) * Delta_t_pulse * omega);
    }
    return result;
}
