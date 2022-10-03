#include "special_curve.h"

#include <cmath>
#include <stdexcept>
#include <iostream>

double Special::xp(double p) {
    return p;
}

double Special::yp(double p) {
    try {
        if (p < pmin || p > pmax) {throw;}
        if (p < -3) {
            return 1.0/2 * 1 / (1 + std::exp(-3*(p+6)));
        } else {
            return 1.0/2 * 1 / (1 + std::exp(3*p));
        }
    } catch (...){
        std::cerr << "Function not defined for x not in [";
        std::cerr << pmin << "," << pmax << "]" << std::endl;
        exit(1);
    }
}

double Special::dxp(double p) {
    return 1;
}

double Special::dyp(double p) {
    try {
        if (p < pmin || p > pmax) {throw;}
        if (p < -3) {
            return 3*std::exp(-3*(p+6))/(2*std::pow(std::exp(-3*(p+6))+1,2));
        } else {
            return -3*std::exp(3*p)/(2*std::pow(std::exp(3*p)+1,2));
        }
    } catch (...){
        std::cerr << "Function not defined for x not in [";
        std::cerr << pmin << "," << pmax << "]" << std::endl;
        exit(1);
    }
}
