#include "curvebase.h"

class Special : public CurveBase {
    public : 
        Special(bool rev=false) : CurveBase(-10, 5, rev) {length = integrate(pmax);} 
    private :
        double xp(double p);
        double yp(double p);
        double dxp(double p);
        double dyp(double p);
};
