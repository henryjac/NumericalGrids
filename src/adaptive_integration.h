#include <cmath>

template <class ASI_fkt>
double simpson_quad(ASI_fkt f, const double a, const double b) {
    return (b-a)/6 * (f(a) + 4*f((a+b)/2) + f(b));
}

template <class ASI_fkt>
double ASI(ASI_fkt f, double a, double b, double tol=1e-8) {
    double c = (a+b)/2;
    double I1 = simpson_quad(f,a,b);
    double I2 = simpson_quad(f,a,c) + simpson_quad(f,c,b);
    double errest = fabs(I1-I2);
    if (errest < 15*tol) {
        return I2;
    }
    return ASI(f,a,(a+b)/2,tol/2) + ASI(f,(a+b)/2,b,tol);
}
