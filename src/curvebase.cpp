#include "curvebase.h"
#include "adaptive_integration.h"

#include <cmath>
#include <iostream>
#include <functional>

CurveBase::CurveBase(double p_min, double p_max, bool _rev) {
    pmin = p_min; pmax = p_max; rev = _rev;
}

// Integrand function for use in Newtons algorithm as
// we need the derivative of the integral.
double CurveBase::integrand(double q) {
    double dxp_q = dxp(q); double dyp_q = dyp(q);
    return std::sqrt(dxp_q*dxp_q + dyp_q*dyp_q);
}

double CurveBase::integrate(double p) {
    // Same as the member functions as I couldnt pass the 
    // member functions to the ASI function
    auto integrand = [&](double q) -> double {
        double dxp_q = dxp(q); double dyp_q = dyp(q);
        return std::sqrt(dxp_q*dxp_q + dyp_q*dyp_q);
    };
    return ASI(integrand,pmin,p); 
}

double CurveBase::find_p(double s) {
    std::function<double(double)> f = [&](double p) -> double {
        return integrate(p) - s*integrate(pmax);
    };
    std::function<double(double)> df = [&](double p) -> double {
        return integrand(p);
    };
    return newton(f, df);
}

Point CurveBase::xy(double s) {
    double p = find_p(s);
    return Point(xp(p),yp(p));
}

// Returns the value of the specific corner, start of the curve if
// true, otherwise the end.
Point CurveBase::get_corner(bool start) {
    if (rev != start) {return xy(1);} else {return xy(0);}
}

bool CurveBase::is_reversed() {
    return rev;
}

void CurveBase::print_corners() {
    Point xy0 = xy(0); Point xy1 = xy(1);
    xy0.print();
    std::cout << " -> ";
    xy1.print('\n');
}

typedef std::function<double(double)> Fnc1D;
double newton(Fnc1D f, Fnc1D df, double x0, double tol, double maxit) {
    int i=0;
    double err=1;
    double x1;
    while (err > tol && i < maxit) {
        x1 = x0 - f(x0)/df(x0);
        err = fabs(x1-x0);
        x0 = x1;
        i++;
    }
    if (err < tol) {return x0;}
    std::cerr << "No convergence" << std::endl;
    exit(1);
}
