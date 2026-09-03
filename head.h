#ifndef __HEAD__
#define __HEAD__

#include </eigen-3.4.0/Eigen/Sparse>
#include <ctime> 

using T = Eigen::Triplet <double>;
using spmat = Eigen::SparseMatrix <double, Eigen::ColMajor>;

/*2dMemory.cpp*/
double** allocate_memory2d(const int M, const int N, const double init_value);
void deallocate_memory2d(double** v);
double*** allocate_memory3d(const int P, const int Q, const int R, double init_value);
void deallocate_memory3d(double*** v, const int P, const int Q);

/*calcWLPcoef.cpp*/
void calcWLPcoef(double* coef, double* array, double** Laguerre);

/*calc_eh.cpp*/
void calc_er_q(double* er_q, double* er_sum, double* Er_q, double* Er_sum);
void calc_hphi_q(double* hphi_q, double* hphi_sum, double* Hphi_q, double* Hphi_sum);

/*calcH.cpp*/
void calc_Hphi(double* Hphi_q, double* Er_q, double* Hphi_sum, double* hphi_sum);

/*composeMatA.cpp*/
void composeMatA (std::vector <T> &t);


/*Fourier.cpp*/
std::complex <double> Fourier(const double* data, double omega, int N);

/*Jr_test.cpp*/
void Jr_test(double* Jq_array, double* Jt_array, double** Laguerre);


/*Jr.cpp*/
double Jr(double t);

/*Laguerre.cpp*/
void makeLaguerreMat(double** Laguerre, double Delta_t, int Nt);
void output_Laguerre(double** Laguerre, double Delta_t, int Nt, std::string filename);

/*PECboundary.cpp*/
void changePECboundary(spmat &A);

/*toTimeDomain.cpp*/
void toTimeDomain(double* timeDomain, double* WLPcoef, double** Laguerre, int Nt);

/*update_beta.cpp*/
void update_beta(double* b, double* Hphi_sum, double* Er_sum, double* hphi_sum, double* er_sum);
void toTimeDomain_Er(int q, double** Er_t, double* Er_q, double** Laguerre);

/*update_sum.cpp*/
void update_sum_r(double* Er_q, double* Er_sum);
void update_sum_phi(double* Hphi_q, double* Hphi_sum);

#endif