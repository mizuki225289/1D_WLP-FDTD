#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <ctime>
#include "const_sq.h"

extern std::string file_dir_name;
///高次のWLP算出(時刻tごとにp次までのWLP φ_p(t) (=v[p]) を計算)

void weighting(const int p, double *v, double *w, double t){
    for(int n = 0; n <= p; n++){
        v[n] *=  std::exp( -w[n] - 0.5*t ); 
    }
}

void wlp(const int p, const double t, double *v){
    /* Weighted Laguerre Polynomial 𝜑ₙ(z) */
    /* v[n] = 𝜑ₙ(t) */

    double *w = new double [p+1];
    for(int i = 0; i <= p; i++){
        w[i] = 0.0;
    }
    
    v[0] = 1.0;

    if ( p > 0 ){
        v[1] = 1.0 - t;

        for(int n = 2; n <= p; n++){
        w[n] = w[n-1];

        /* exp(10)を超える毎に exp(-10)倍する。何倍したか w に記憶しておく */
        if ( std::log(std::abs(v[n-1])) > 10.0 ){
        v[n-1] *= std::exp(-10.0);
        v[n-2] *= std::exp(-10.0);
        w[n] -= 10.0;
        w[n-1] -= 10.0;
        w[n-2] -= 10.0;
        }

        v[n] = (2.0*n - 1.0 - t)/n * v[n-1] - (n - 1.0)/n * v[n-2];
        }
    }

    /* exp(-10 * N)倍を戻すのと重み exp(-t/2) を同時に戻す */
    weighting(p, v, w, t);
    delete [] w;

}

void makeLaguerreMat(double** Laguerre, double Delta_t, int Nt) {
    double v[q_MAX+1];
    for(int i=0; i < Nt; i++) {
        double t = i * Delta_t;
        wlp(q_MAX, s * t, v);
        for(int q=0; q < q_MAX; q++) {
            Laguerre[q][i] = v[q];
        }
    }
}

void output_Laguerre(double** Laguerre, double Delta_t, int Nt, std::string filename) {
    std::ofstream ofs_Laguerre(file_dir_name + filename + ".dat");
    for(int t=0; t < Nt; t++) {
        ofs_Laguerre << t * Delta_t << " ";
        for(int p=0; p < q_MAX; p++) {
            ofs_Laguerre << Laguerre[p][t] << " ";
        }
        ofs_Laguerre << std::endl;
    }
    ofs_Laguerre.close();
}