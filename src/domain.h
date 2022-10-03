#include "curvebase.h"

#include <cmath>
#include <memory>

enum StretchDir {
    STRETCH_X, STRETCH_Y
};

class Domain {
    static double phi0(double s) {return 1-s;}
    static double phi1(double s) {return s;}
    private: 
        std::shared_ptr<CurveBase> borders[4];
        int n,m;
        std::unique_ptr<double[]> x,y;
        std::function<double(double)> stretch_x, stretch_y;
        bool check_consistency(); // Checks that all the curves end where the next starts.
    public:
        Domain(std::shared_ptr<CurveBase>, std::shared_ptr<CurveBase>, 
                std::shared_ptr<CurveBase>, std::shared_ptr<CurveBase>);
        Domain(const Domain&);
        Domain& operator=(const Domain&);
        ~Domain() {};
        Point operator()(int, int);
        
        void set_stretching(std::function<double(double)>, StretchDir);
        void generate_grid(int n_, int m_);
        void print_corners();

        int xsize();
        int ysize();
        double getx(int, int);
        double gety(int, int);

        void save_boundary(const char*, int precision=50);
        void save_grid(const char*);
};
