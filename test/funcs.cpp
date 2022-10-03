#include "../src/special_curve.h"
#include "../src/straight_line.h"
#include "../src/grid_functions.h"

#include <iostream>
#include <array>
#include <memory>
#include <cmath>

void calc_and_save(bool w_stretching) {
    std::shared_ptr<Special> sp = std::make_shared<Special>();
    std::shared_ptr<StraightLine> s1 = std::make_shared<StraightLine>(0,1,5,0,0,3);
    std::shared_ptr<StraightLine> s2 = std::make_shared<StraightLine>(1,0,0,3,-10,5,true);
    std::shared_ptr<StraightLine> s3 = std::make_shared<StraightLine>(0,1,-10,0,0,3,true);

    std::shared_ptr<Domain> d = std::make_shared<Domain>(sp,s1,s2,s3); // Normal

    if (w_stretching) {
        double delta(3);
        std::function<double(double)> stretch = [delta](double sigma) -> double {
            return 1 + std::tanh(delta*(sigma-1)) / std::tanh(delta);
        };
        d->set_stretching(stretch, STRETCH_Y);
    }

    d->generate_grid(35,30);

    GFkt U(d);
    GFkt U_analyticpdx(d);
    GFkt U_analyticpdy(d);
    GFkt U_analyticlaplace(d);

    // Analytic functions
    auto fnc = [](Point P) -> double {
        double t = P.x()/10;
        return std::sin(t*t)*std::cos(t)+P.y();
    };
    auto fnc_pdx = [](Point P) -> double {
        double t = P.x()/10;
        return std::cos(t*t)*t*std::cos(t)/5 - std::sin(t*t)*sin(t)/10;
    };
    auto fnc_pdy = [](Point P) -> double {
        return 1;
    };
    auto fnc_laplace = [](Point P) -> double {
        double t =  P.x()/10;
        return (-100*t*std::sin(t)*std::cos(t*t) - std::cos(t)*((100*t*t+25)*std::sin(t*t) -
             50*std::cos(t*t)) ) / 2500;
    };

    // Fill with analytic values
    U.fill_matrix(fnc);
    U_analyticpdx.fill_matrix(fnc_pdx);
    U_analyticpdy.fill_matrix(fnc_pdy);
    U_analyticlaplace.fill_matrix(fnc_laplace);

    std::array<std::shared_ptr<GFkt>,2> Updxy = U.pd();
    GFkt Updx = *Updxy[0];
    GFkt Updy = *Updxy[1];
    GFkt Ulaplace = U.laplace();

    // Calculate error norms
    double pdx_err = (Updx - U_analyticpdx).norm_1();
    double pdy_err = (Updy - U_analyticpdy).norm_1();
    double laplace_err = (Ulaplace - U_analyticlaplace).norm_1();

    // Print errors
    if (!w_stretching) {
        std::cout << "Uniform grid" << std::endl;
    } else {
        std::cout << "Nonuniform grid" << std::endl;
    }
    std::cout << "------------------------" << std::endl;
    std::cout << "pdx error:\t" << pdx_err << std::endl;
    std::cout << "pdy error:\t" << pdy_err << std::endl;
    std::cout << "laplace error:\t" << laplace_err << std::endl;
    std::cout << "------------------------" << std::endl;

    // Save binary files
    std::string grid = "bin/data/grid";
    std::string boundary = "bin/data/boundary";
    std::string function = "bin/data/function";
    std::string analyticpdx = "bin/data/analyticpdx";
    std::string analyticpdy = "bin/data/analyticpdy";
    std::string analyticlaplace = "bin/data/analyticlaplace";
    std::string pdx = "bin/data/pdx";
    std::string pdy = "bin/data/pdy";
    std::string laplace = "bin/data/laplace";

    /* std::array<std::string,9> strings = { */
    std::string *strings[9] = {
        &grid,&boundary,&function,
        &analyticpdx,&analyticpdy,&analyticlaplace,
        &pdx,&pdy,&laplace
    };

    std::string add;
    if (!w_stretching) {
        add = ".bin";
    } else {
        add ="_stretch.bin";
    }
    for (int i=0; i<9; i++) {
        *strings[i] += add;
    }

    U.save(function.c_str(), grid.c_str(), boundary.c_str());

    U_analyticpdx.save(analyticpdx.c_str());
    U_analyticpdy.save(analyticpdy.c_str());
    U_analyticlaplace.save(analyticlaplace.c_str());

    Updx.save(pdx.c_str());
    Updy.save(pdy.c_str());
    Ulaplace.save(laplace.c_str());
}

int main() {
    calc_and_save(false);
    calc_and_save(true);
}
