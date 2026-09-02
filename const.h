#ifndef __CONST__
#define __CONST__

#define _USE_MATH_DEFINES
#include <math.h>
#include <string>

/*定数*/
constexpr double C0 {2.99792458e8}; /*真空中の光速*/
constexpr double MU0 {4.0 * M_PI * 1.0e-7}; /*真空中の透磁率*/
constexpr double EPS0 {1.0 / MU0 / C0 / C0}; /*真空中の誘電率*/
constexpr double e0 {1.60e-19}; /*電気素量*/
constexpr double me {9.10e-31}; /*電子の質量*/
constexpr double EARTH_R {6378.1e3}; //地球の赤道半径

#endif