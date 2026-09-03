#include <vector>
#include <fstream>
#include <iostream>
#include "head.h"
#include "const.h"
#include "const_sq.h"
#include "const_region.h"
#include "r_theta.h"
#include "PMLparameter.h"

int count = 0;
constexpr double a0 {2.0 / EPS0 / s};
constexpr double b0 {2.0 / MU0 / s};
constexpr double xi0 {a0 * b0}; // 4.0 / EPS0 / MU0 / s / s

/*中心軸　j=0*/
void centerAxis(std::vector <T> &t) {
    int idx = 0;

    int Center = 0;
    int Right = 1;

    t.push_back( T(idx, Center, 1.0 + 4.0 * xi0 / EARTH_R / EARTH_R / Delta_THETA / Delta_THETA) );
    t.push_back( T(idx,  Right,     - 4.0 * xi0 / EARTH_R / EARTH_R / Delta_THETA / Delta_THETA) );
    count++;
}

void center(std::vector <T> &t) {
    for(int j=1; j < N_THETA - L; j++) {
        int idx = j;

        int Left = j - 1;
        int Center = j;
        int Right = j + 1;

        double C1 = xi0 / EARTH_R / EARTH_R / sin_j(j) / Delta_THETA / Delta_THETA;

        t.push_back( T(idx, Center, 1.0 + C1 * (sin_j(j+0.5) + sin_j(j-0.5))) );
        t.push_back( T(idx,  Right,     - C1 * sin_j(j+0.5) ) );
        t.push_back( T(idx,   Left,     - C1 * sin_j(j-0.5)) );
        count++;
    }
}

void PML(std::vector <T> &t) {
    for(int j = N_THETA - L; j < N_THETA; j++) {
        int idx = j;

        int Left = j - 1;
        int Center = j;
        int Right = j + 1;

        double C1 = xi_theta(j) * xi0 / EARTH_R / EARTH_R / sin_j(j) / Delta_THETA / Delta_THETA;

        t.push_back( T(idx, Center, 1.0 + C1 * (sin_j(j+0.5) * xi_theta(j+0.5) + sin_j(j-0.5) * xi_theta(j-0.5)) ) );
        t.push_back( T(idx,  Right,     - C1 * sin_j(j+0.5) * xi_theta(j+0.5) ) );
        t.push_back( T(idx,   Left,     - C1 * sin_j(j-0.5) * xi_theta(j-0.5)) );
        count++;
    }
}

/*右端　j = N_THETA*/
void right_edge(std::vector <T> &t) {
    int idx = N_THETA;

    int Left = N_THETA - 1;
    int Center = N_THETA;

    t.push_back( T(idx,   Left, 0.0) );
    t.push_back( T(idx, Center, 1.0) );
    count++;
}

void composeMatA(std::vector <T> &t) {
    centerAxis(t);
    center(t);
    PML(t);
    right_edge(t);

    std::cout << "count = " << count << ", Matsize = " << Matsize << std::endl;
}
