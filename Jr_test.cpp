#include <iostream>
#include <fstream>
#include <string>
#include <ctime>

#include "head.h"
#include "const_sq.h"
#include "const_time.h"

extern std::string file_dir_name;

void Jr_test(double* Jq_array, double* Jt_array, double** Laguerre) {
    /*Iq_arrayのテスト1　係数をそのままファイルに書き出す*/
    std::ofstream ofs_Jq_test1(file_dir_name + "Jq_array_test.dat");
    for(int i=0; i < q_MAX; i++) {
        ofs_Jq_test1 << i << " " << Jq_array[i] << std::endl;
    }
    ofs_Jq_test1.close();

    /*Iq_arrayのテスト2　復元した結果(It_calc)と元の式で計算した結果(It_array)を比較*/
    double* Jt_calc = new double [Nt_pulse];
    toTimeDomain(Jt_calc, Jq_array, Laguerre, Nt_pulse);

    std::ofstream ofs_Jq_test2(file_dir_name + "Jr_compare.dat");
    ofs_Jq_test2 << "# time origin calc_fromWLP" << std::endl;
    for(int i=0; i < Nt_pulse; i++) {
        ofs_Jq_test2 << i * Delta_t_pulse << " " << Jt_array[i] << " " << Jt_calc[i] << std::endl;
    }
    ofs_Jq_test2.close();
    delete [] Jt_calc;
    std::cout << "Iq test OK" << std::endl;
}