#ifndef __CONST_OPTION__
#define __CONST_OPTION__

struct SimulationOptions {
    /*電流印加点*/
    bool isJr_j0all = true; /*中心軸全体に電流を流すときtrue*/
    int input_i = 0; /*input_i, input_j は isJr_j0all = false のときに利用*/
    int input_j = 0;

    /*プラズマ領域*/
    double omegap_uni = 0.0; /*↑trueのときのプラズマ角周波数[rad/s]*/
    double nu_uni = 0.0; /*↑trueのときの衝突周波数[Hz]*/

    /*観測点設定*/
    int obs_i = 25;
    int obs_j = 100;
};

extern SimulationOptions opt;

#endif