#ifndef CURVEBASE
#define CURVEBASE

#include "point.h"

#include <functional>

class CurveBase {
    protected :
        double pmin;
        double pmax;
        double length;
        bool rev;

        virtual double xp(double p) = 0;
        virtual double yp(double p) = 0;
        virtual double dxp(double p) = 0;
        virtual double dyp(double p) = 0;
        double integrand(double q);
        double integrate(double p);
        double find_p(double s);
    public :
        CurveBase(double p_min=0, double p_max=1, bool rev=false);
        virtual ~CurveBase() {};
        Point xy(double s); // To get x/y do .first/.second
        Point get_corner(bool start);
        bool is_reversed();
        void print_corners();
};

typedef std::function<double(double)> Fnc1D;
double newton(Fnc1D f, Fnc1D df, double x0=0, double tol=1e-12, double maxit=1000);
#endif
