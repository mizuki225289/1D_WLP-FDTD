#ifndef __CONST_TIME__
#define __CONST_TIME__

#include "const_region.h"

constexpr double t_MAX_pulse {1.0e-3}; /*印加電流が十分落ち着くまでの時間*/
constexpr int Nt_pulse {10000};
constexpr double Delta_t_pulse {t_MAX_pulse / Nt_pulse};

constexpr double t_MAX_output {2.0e-3}; /*この時刻まで出力する 波源からθ方向の端点まで行くまで*/
constexpr int Nt_output {500};
constexpr double Delta_t_output {t_MAX_output / Nt_output};

#endif