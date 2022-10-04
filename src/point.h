#ifndef POINT 
#define POINT 
#include <iostream>

#include <cmath>

class Point {
    private:
        double x_;
        double y_;
    public:
        Point() : x_(0), y_(0) {};
        Point(double x, double y) : x_(x), y_(y) {}
        friend std::ostream& operator<<(std::ostream& os, const Point& P);
        Point operator+(Point p) {return Point(x_+p.x_, y_+p.y_);}
        Point operator-(Point p) {return Point(x_-p.x_, y_-p.y_);}
        double x() const {return x_;}
        double y() const {return y_;}
        void print(char end='\0');
};

Point operator*(Point a, double b);
Point operator*(double a, Point b);
bool eps_equal(Point,Point,double eps=1e-10);

#endif
