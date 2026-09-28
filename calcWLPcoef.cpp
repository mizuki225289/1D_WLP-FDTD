#include "const_time.h"
#include "const_sq.h"

double integral(double* p) {
    /*シンプソンの公式を用いた積分*/
    double result = 0.0;
    double dt = Delta_t_pulse;
    for(int i=0; i < Nt_pulse / 2; i++) {
        result += s * (dt / 3.0) * (p[2*i] + 4.0 * p[2*i+1] + p[2*i+2]);
    }

    return result;
}

void calcWLPcoef(double* coef, double* array, double** Laguerre) {
    double *tmp = new double [Nt_pulse + 1];
    for(int q=0; q < q_MAX; q++) {
        for(int i=0; i < Nt_pulse + 1; i++) {
            tmp[i] = Laguerre[q][i] * array[i];
        }
        coef[q] = integral(tmp);
    }
    delete[] tmp;
}