#include "straight_line.h"

double StraightLine::xp(double p) {
    return a*p+c;
}

double StraightLine::yp(double p) {
    return b*p+d;;
}

double StraightLine::dxp(double p) {
    return a;
}

double StraightLine::dyp(double p) {
    return b;
}
