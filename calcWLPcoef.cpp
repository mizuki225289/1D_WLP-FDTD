#include "const_time.h"
#include "const_sq.h"

double integral(double* p, int N) {
    /*シンプソンの公式を用いた積分*/
    double result = 0.0;
    double dt = t_MAX_pulse / N;
    for(int i=0; i < N/2; i++) {
        result += s * (dt / 3.0) * (p[2*i] + 4.0 * p[2*i+1] + p[2*i+2]);
    }

    return result;
}

void calcWLPcoef(double* coef, double* array, double** Laguerre) {
    for(int q=0; q < q_MAX; q++) {
        double *tmp = new double [Nt_pulse];
        for(int i=0; i < Nt_pulse; i++) {
            tmp[i] = Laguerre[q][i] * array[i];
        }
        coef[q] = integral(tmp, Nt_pulse);
    }
}