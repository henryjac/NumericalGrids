#include "../src/straight_line.h"
#include "../src/special_curve.h"
#include "../src/domain.h"

#include <iostream>
#include <iomanip>
#include <memory>

void test_curves();
void test_domain();

int main() {
    test_curves();
    test_domain();
}

void test_curves() {
    Special B;
    double between[11] = {0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1};
    Point xy;
    std::cout << "Special curve arc length parametrization" << std::endl;
    for (int i=0; i<11; i++) {
        xy = B.xy(between[i]);
        std::cout << std::setprecision(5) << std::fixed;
        std::cout << xy << std::endl;
    }
}

void test_domain() {
    std::shared_ptr<Special> sp = std::make_shared<Special>();
    std::shared_ptr<StraightLine> s1 = std::make_shared<StraightLine>(0,1,5,0,0,3);
    std::shared_ptr<StraightLine> s2 = std::make_shared<StraightLine>(1,0,0,3,-10,5,true);
    std::shared_ptr<StraightLine> s3 = std::make_shared<StraightLine>(0,1,-10,0,0,3,true);

    Domain d(sp,s1,s2,s3);
    Domain D(d);

    int n=50,m=20;

    d.generate_grid(n,m);

    char saveboundary[] = "bin/data/boundary.bin";
    char savegrid[] = "bin/data/grid.bin";
    d.save_boundary(saveboundary, 50);
    d.save_grid(savegrid);

    double delta(3);
    std::function<double(double)> stretch = [delta](double sigma) -> double {
        return 1 + std::tanh(delta*(sigma-1)) / std::tanh(delta);
    };
    D.set_stretching(stretch, STRETCH_Y);
    D.generate_grid(n,m);
    char savegrid_stretched[] = "bin/data/grid_stretched.bin";
    D.save_grid(savegrid_stretched);
}
