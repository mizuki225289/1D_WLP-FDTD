#ifndef __CONST_REGION__
#define __CONST_REGION__

#define _USE_MATH_DEFINES
#include <math.h>
#include "const.h"

/*シミュレーション領域について*/
constexpr double arc_length {300e3}; /*θ方向の孤の長さ（地表面）*/
constexpr double THETA_MAX {arc_length / EARTH_R};
constexpr int N_THETA {300};
constexpr double Delta_THETA {THETA_MAX / N_THETA};
constexpr double Delta_arc {arc_length / N_THETA};
constexpr int Matsize {N_THETA + 1};

/*PMLについて*/
constexpr int L {10}; //PMLのレイヤー数
constexpr int K {3}; /*PML領域のちょっと内側から補助量を計算しておく*/

#endif