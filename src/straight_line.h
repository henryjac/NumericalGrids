#include "curvebase.h"

class StraightLine : public CurveBase {
    public : 
        StraightLine(
                double aa, double bb, double cc, double dd, 
                double p_min=0, double p_max=1, bool rev=false)
                : CurveBase(p_min, p_max, rev) {
            a=aa; b=bb; c=cc; d=dd;
        }
        StraightLine() : a(), b(), c(), d() {};
    private :
       double a, b, c, d; 
       // [x,y] = [a,b]t + [c,d]
       // As opposed to (ax + by + c = 0) as that is not trivial to parametrize so would need some conditional statements
       double xp(double p);
       double yp(double p);
       double dxp(double p);
       double dyp(double p);
};
