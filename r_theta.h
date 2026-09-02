#ifndef __R_THETA__
#define __R_THETA__

#include <cmath>
#include "const.h"
#include "const_region.h"

template <typename ValueType>
double theta_j(ValueType j) {
    return Delta_THETA * j;
}

template <typename ValueType>
double sin_j(ValueType j) {
    return sin(Delta_THETA * j);
}

#endif