#include "const_region.h"
#include "head.h"

void changePECboundary(spmat &A) {
    int idx = N_THETA;
    int Left = N_THETA - 1;

    A.coeffRef(idx, idx) = 1.0;
    A.coeffRef(idx,   Left) = 0.0;
    A.coeffRef(Left, idx) = 0.0;
}