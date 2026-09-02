// #include </eigen-3.4.0/Eigen/SparseLU>
#include </eigen-3.4.0/Eigen/IterativeLinearSolvers>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <ctime>
#include <filesystem>
#include <chrono>

#include "head.h"
#include "const.h"
#include "const_time.h"
#include "const_sq.h"
#include "const_region.h"
#include "option.h"
#include "r_theta.h"
#include "PMLparameter.h"

using Solver = Eigen::SparseLU < spmat, Eigen::COLAMDOrdering <int> >; /*LU分解*/
std::string file_dir_name;

SimulationOptions opt;

int main(void) {
    time_t t;
    t = time(NULL);

    auto now = std::chrono::system_clock::now();
    auto now_c = std::chrono::system_clock::to_time_t(now);
    std::stringstream ss;
    ss << std::put_time(localtime(&now_c), "%Y%m%d_%H%M");
    std::cout << ss.str() << " " << std::ctime(&t);

    // フォルダを作成
    std::filesystem::create_directory("./result_" + ss.str());
    file_dir_name = "./result_" + ss.str() + "/";

    /*所要時間計測用*/
    std::chrono::system_clock::time_point start, end;
    double elapsed;

    /*Laguerre作成*/
    std::cout << "Make Laguerre" << std::endl;
    start = std::chrono::system_clock::now();
        double** Laguerre_pulse = allocate_memory2d(q_MAX, Nt_pulse, 0.0); //電流パルス用の時間幅が細かいWLP
        makeLaguerreMat(Laguerre_pulse, Delta_t_pulse, Nt_pulse);
        // output_Laguerre(Laguerre_pulse, Delta_t_pulse, Nt_pulse, "Laguerre_pulse");

        double** Laguerre_output = allocate_memory2d(q_MAX, Nt_output, 0.0); //時間波形出力用の時間幅が大きいWLP
        makeLaguerreMat(Laguerre_output, Delta_t_output, Nt_output);
        // output_Laguerre(Laguerre_output, Delta_t_output, Nt_output, "Laguerre_output");
    end = std::chrono::system_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::milliseconds> (end - start).count();
    std::cout << "Laguerre OK  elapsed = " << elapsed << " msec." << std::endl;

    /*電流印加点設定*/
    int input_index = 0;

    /*Jt_array作成*/
    double* Jt_array = new double [Nt_pulse];
    for(int i=0; i < Nt_pulse; i++) {
        Jt_array[i] = Jr(i * Delta_t_pulse);
    }
    double* Jq_array = new double [q_MAX];
    calcWLPcoef(Jq_array, Jt_array, Laguerre_pulse);
    // std::cout << "bbb" << std::endl;

    /*Irの結果を確認する（WLPの展開係数、WLPからの復元結果と式から計算した値の比較）*/
    Jr_test(Jq_array, Jt_array, Laguerre_pulse);

    /*Ideltal テスト*/
    std::ofstream ofs_Ideltal(file_dir_name + "Ideltal(t).dat");
    for(int i=0; i < Nt_pulse; i++) {
        ofs_Ideltal << i * Delta_t_pulse << " " << Jt_array[i] << std::endl;
    }
    ofs_Ideltal.close();

    /*Irのフーリエ変換*/
    std::ofstream ofs4(file_dir_name + "Iomega.dat");
    for(int i=0; i < 600; i++) { /*10^3 ~ 10^6*/
        double omega = pow(10.0, 3.0 + i / 200.0);
        ofs4 << omega << " " << abs(Fourier(Jt_array, omega, Nt_pulse)) << std::endl;
    }
    ofs4.close();

    /* Tripletの用意 */
    std::vector <T> CoefMatrix;
    CoefMatrix.reserve( 3 * Matsize ); /* 予め、記憶領域を確保しておく */
    composeMatA( CoefMatrix );
    // std::cout << "cc" << std::endl;

    /* 疎行列の用意 */
    std::cout << "Make SparseMat" << std::endl;
    start = std::chrono::system_clock::now();
        spmat A(Matsize, Matsize);
        A.setFromTriplets( CoefMatrix.begin(), CoefMatrix.end() );
        changePECboundary(A);
        A.makeCompressed();
    end = std::chrono::system_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::milliseconds> (end - start).count();
    std::cout << "SparseMat OK  elapsed = " << elapsed << " msec." << std::endl;
    std::cout << "ccc" << std::endl;

    /*LU分解利用　ここから*/
        /* LU分解 */
        std::cout << "LU decomposition" << std::endl;
        Solver sol;
        start = std::chrono::system_clock::now();
            sol.analyzePattern( A );
            sol.factorize( A );
        end = std::chrono::system_clock::now();
        elapsed = std::chrono::duration_cast<std::chrono::seconds> (end - start).count();
        std::cout << "factorize OK  elapsed = " << elapsed << " sec. = " << elapsed / 60.0 << "min." << std::endl; /*57万×57万で約15分*/

        /*LU分解成功判定*/
        if(sol.info() == Eigen::Success) {
            std::cout << "Success!!" << std::endl;
            // exit(0);
        } else {
            std::cout << sol.lastErrorMessage() << std::endl;
            exit(1);
        }
    /*LU分解　ここまで*/

    /*配列確保*/
    double* Er_q = new double [N_THETA + 1];
    double* Er_sum = new double [N_THETA + 1];
    double* Hphi_q = new double [N_THETA];
    double* Hphi_sum = new double [N_THETA];

    double* er_q = new double [N_THETA + 1];
    double* er_sum = new double [N_THETA + 1];
    double* hphi_q = new double [N_THETA];
    double* hphi_sum = new double [N_THETA];

    /*初期化*/
    for(int j=0; j < N_THETA+1; j++) {
        Er_sum[j] = 0;
        er_sum[j] = 0;
        Hphi_sum[j] = 0;
        hphi_sum[j] = 0;
    }

    /*Hphiのアニメーション作成用*/
    double** Er_t   = allocate_memory2d(N_THETA + 1, Nt_output, 0.0);

    /*右辺ベクトルb, 解のベクトルx*/
    double* beta = new double [Matsize]; /*b^q = J^q + tilde{beta}*/
    Eigen::VectorXd b(Matsize), x(Matsize);
    for(int i=0; i < Matsize; i++) {
        beta[i] = 0.0;
        b(i) = 0.0;
        x(i) = 0.0;
    }

    /*観測点用意*/
    int obs_j = 100;
    double* obs_Eq = new double [q_MAX];

    std::cout << "main" << std::endl;
    start = std::chrono::system_clock::now();
    std::chrono::system_clock::time_point mid = std::chrono::system_clock::now();

    /*メインループ*/
    for(int q=0; q < q_MAX; q++) {
        if(q%50 == 0) {
            end = std::chrono::system_clock::now();
            elapsed = std::chrono::duration_cast<std::chrono::milliseconds> (end - mid).count();
            std::cout << "q = " << q << " / " << q_MAX << ", elapsed = " << elapsed << " msec." << std::endl;
            mid = std::chrono::system_clock::now();
        }

        /*右辺更新*/
        for(int i=0; i < Matsize; i++) {
            b[i] = beta[i];
        }
        b[input_index] += (-2.0 / EPS0 / s) * Jq_array[q];

        /*方程式解くbyLU*/
        x = sol.solve(b);

        /*エラー検出*/
        if (sol.info() != Eigen::Success) {
            std::cout << "solver failed at q = " << q << std::endl;
            exit(1);
        }
        if (!x.allFinite()) {
            std::cout << "NaN or Inf detected at q = " << q << std::endl;
            exit(1);
        }

        /*x から Er_q, Etheta_q作成*/
        for(int j=0; j < N_THETA + 1; j++) {
            Er_q[j] = x[j];
        }

        /*Hphi_q, H_phir_q, H_phitheta_q 計算*/
        calc_Hphi(Hphi_q, Er_q, Hphi_sum, hphi_sum);

        /*e_theta*/
        calc_er_q(er_q, er_sum, Er_q, Er_sum);
        calc_hphi_q(hphi_q, hphi_sum, Hphi_q, Hphi_sum);

        /*sum 更新*/
        update_sum_r(Er_q, Er_sum);
        update_sum_r(er_q, er_sum);
        update_sum_phi(Hphi_q, Hphi_sum);
        update_sum_phi(hphi_q, hphi_sum);
        // std::cout << "iii" << std::endl;

        /*右辺更新*/
        update_beta(beta, Hphi_sum, Er_sum, hphi_sum, er_sum);
        // std::cout << "jjj" << std::endl;

        /*時間波形計算（逐次更新）*/
        toTimeDomain_Er(q, Er_t, Er_q, Laguerre_output);
        // std::cout << "kkk" << std::endl;

        /*観測点の展開係数保存*/
        obs_Eq[q] = Er_q[obs_j];

        // std::cout << "Er^q = " << Er_q[295] << std::endl;
        // std::cout << "Hphi^q = " << Hphi_q[295] << std::endl;
        // std::cout << "Er_sum = " << Er_sum[295] << std::endl;
        // std::cout << "Hphi_sum = " << Hphi_sum[295] << std::endl;
        // std::cout << "er^q = " << er_q[295] << std::endl;
        // std::cout << "hphi^q = " << hphi_q[295] << std::endl;
        // std::cout << "er_sum = " << er_sum[295] << std::endl;
        // std::cout << "hphi_sum = " << hphi_sum[295] << std::endl;
    }
    end = std::chrono::system_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::milliseconds> (end - start).count();
    std::cout << "Main loop comp.  elapsed = " << elapsed << " msec. = " << elapsed / 60.0 << "sec." << std::endl;
    std::cout << elapsed / q_MAX << " msec. per once loop." << std::endl;

    /*アニメーション作成用*/
    std::filesystem::create_directory(file_dir_name + "Er_t");
    for(int t=0; t < Nt_output; t++) {
        if(t % 100 == 0) {
            std::cout << t << " / " << Nt_output << std::endl;
        }
        std::ofstream ofs_animation(file_dir_name + "Er_t/t_" + std::to_string(t) + ".dat");
        for(int j=0; j < N_THETA + 1; j++) {
            ofs_animation << Delta_arc * j * 1.0e-3 << " " << Er_t[j][t] << std::endl;
        }
    }
    std::cout << "complete!" << std::endl;
    deallocate_memory2d(Er_t);

    /*観測点の時間波形計算*/
    double* obs_Et = new double [Nt_output];
    toTimeDomain(obs_Et, obs_Eq, Laguerre_output, Nt_output);
    std::ofstream ofs_obs(file_dir_name + "obs_Et.dat");
    for(int i=0; i < Nt_output; i++) {
        ofs_obs << i * Delta_t_output << " " << obs_Et[i] << std::endl;
    }
    delete[] obs_Eq;
    delete[] obs_Et;

    delete[] Er_q;
    delete[] Er_sum;
    delete[] Hphi_q;
    delete[] Hphi_sum;
    delete[] er_q;
    delete[] er_sum;
    delete[] hphi_q;
    delete[] hphi_sum;
    deallocate_memory2d(Laguerre_output);
    deallocate_memory2d(Laguerre_pulse);

    delete [] Jt_array;
    delete [] Jq_array;
    delete [] beta;

    return 0;
}