#include "point.h"

void Point::print(char end) {
    std::cout << "(" << x_ << "," << y_ << ")" << end;
}
 Point operator*(Point a, double b) {
    return Point(a.x()*b, a.y()*b);
}

Point operator*(double a, Point b) {
    return Point(b.x()*a,b.y()*a);
}

bool eps_equal(Point a, Point b, double eps) {
    if (
        fabs(a.x() - b.x()) < eps 
        && 
        fabs(a.y() - b.y()) < eps
    ) return true;
    else return false;
}

std::ostream& operator<<(std::ostream& os, const Point& P) {
    return os << "(" << P.x() << "," << P.y() << ")";
}
