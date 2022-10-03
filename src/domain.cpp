#include "domain.h"

#include <iostream>
#include <cmath>
#include <cstdio>

bool Domain::check_consistency() {
    Point first,second;
    for (int i=0; i<4; i++) {
        first = borders[i]->get_corner(true); second = borders[(i+1)%4]->get_corner(false);
        if (!(eps_equal(first, second, 1e-5))) {
            return false; 
        }
    }
    return true;
}

Domain::Domain(std::shared_ptr<CurveBase> c0, std::shared_ptr<CurveBase> c1, 
        std::shared_ptr<CurveBase> c2, std::shared_ptr<CurveBase> c3) {
    n = 0; m = 0; x = nullptr; y = nullptr;
    borders[0] = c0; borders[1] = c1; 
    borders[2] = c2; borders[3] = c3; 

    stretch_x = [](double t) -> double {return t;};
    stretch_y = [](double t) -> double {return t;};

    if (!check_consistency()) {
        std::cout << "Boundary curves does not form a closed surface." << std::endl;
        for (int i=0; i<4; i++) {
            borders[i] = nullptr;
        }
        exit(1);
    }
}

Domain::Domain(const Domain& D) : n(D.n), m(D.m), x(nullptr), y(nullptr), 
        stretch_x(D.stretch_x), stretch_y(D.stretch_y) {
    for (int i=0; i<4; i++) {
        borders[i] = D.borders[i];
    }
    if (m > 0) {
        x = std::make_unique<double[]>(n*m);
        y = std::make_unique<double[]>(n*m);
        for (int i=0; i<n*m; i++) {
            x[i] = D.x[i];
            y[i] = D.y[i];
        }
    }
}

Domain& Domain::operator=(const Domain& D) {
    if (this != &D) {
        if (n*m != D.n*D.m) {
            if (x != nullptr || y != nullptr) {
                x = nullptr; y = nullptr;
            } 
            if (D.x != nullptr || D.y != nullptr) {
                x = std::make_unique<double[]>(D.n*D.m);
                y = std::make_unique<double[]>(D.n*D.m);
            }
        } 
        n = D.n;
        m = D.m;
        for (int i=0; i<n*m; i++) {
            x[i] = D.x[i];
            y[i] = D.y[i];
        }
    }
    return *this;
}

Point Domain::operator()(int i, int j) {
    if (i < 0 || i >= n || j < 0 || j >= m) {
        std::cerr << "Index out of range." << std::endl;
        exit(1);
    }
    return Point(x[i*m+j],y[i*m+j]);
}

void Domain::set_stretching(std::function<double(double)> s, StretchDir dir) {
    switch(dir) {
        case(STRETCH_X): {
            stretch_x = s;
            break;
        }
        case(STRETCH_Y): {
            stretch_y = s;
            break;
        }
    }
}

void Domain::generate_grid(int n_, int m_) {
    if ( n_ < 0 || m_ < 0 ) {
        std::cerr << "Grid generation with negative values not allowed." << std::endl;
        exit(1);
    }

    n = n_; m = m_;
    x.reset(new double[n*m]);
    y.reset(new double[n*m]);

    Point edges;
    Point corners;
    Point xy;
    
    Point *borders0 = new Point[n];
    Point *borders1 = new Point[m];
    Point *borders2 = new Point[n];
    Point *borders3 = new Point[m];

    double *xis = new double[n];
    double *etas = new double[m];

    // Calculate the xy(eta) and xy(xi) to reduce time spent calculating them
    // in the nested for loop
    for (int i=0; i<n; i++) {
        xis[i] = stretch_x((double)i/(n-1)); // -1 to get to 1
        borders0[i] = (*borders[0]).xy(xis[i]);
        borders2[i] = (*borders[2]).xy(xis[i]);
    }
    for (int j=0; j<m; j++) {
        etas[j] = stretch_y((double)j/(m-1));
        /* if (stretching) { etas[j] = stretch(etas[j]); } */
        borders1[j] = (*borders[1]).xy(etas[j]);
        borders3[j] = (*borders[3]).xy(etas[j]);
    }

    // Use algebraic grid generation to generate grid using interpolation
    for (int i=0; i<n; i++) {
        for (int j=0; j<m; j++) {
            edges = phi0(xis[i])*borders3[j] 
                + phi1(xis[i])*borders1[j]
                + phi0(etas[j])*borders0[i] 
                + phi1(etas[j])*borders2[i];
            corners = -phi0(xis[i])*phi0(etas[j])*borders0[0] 
                - phi0(xis[i])*phi1(etas[j])*borders2[0]
                - phi0(etas[j])*phi1(xis[i])*borders0[n-1] 
                - phi1(etas[j])*phi1(xis[i])*borders2[n-1];
            xy = edges + corners;
            x[i*m+j] = xy.x();
            y[i*m+j] = xy.y();
        }
    }

    delete[] borders0; delete[] borders1;
    delete[] borders2; delete[] borders3;
    delete[] xis; delete[] etas;
}

void Domain::print_corners() {
    for (int i=0; i<4; i++) {
        borders[i]->print_corners();
    }
}

int Domain::xsize() {
    return n;
}

int Domain::ysize() {
    return m;
}

double Domain::getx(int i, int j) {
    return x[i*m+j];
}

double Domain::gety(int i, int j) {
    return y[i*m+j];
}

void Domain::save_boundary(const char* file, int precision) {
    FILE* pfile = fopen(file, "wb");
    // Write the precision so julia can determine it explicitly
    fwrite(&precision, sizeof(int), 1, pfile); 

    Point xy;
    double x,y;
    double s;
    for (int i=0; i<4; i++) {
        for (double j=0; j<=precision; j++) {
            if (borders[i]->is_reversed()) {s = 1 - j/precision;}
            else {s = j/precision;}
            xy = borders[i]->xy(s);
            x = xy.x(); y = xy.y();
            fwrite(&x, sizeof(double), 1, pfile);
            fwrite(&y, sizeof(double), 1, pfile);
        }
    }
    fclose(pfile);
}

void Domain::save_grid(const char* file) {
    FILE* pfile = fopen(file, "wb");
    // Write the shape so julia can determine it explicitly
    fwrite(&n, sizeof(int), 1, pfile);
    fwrite(&m, sizeof(int), 1, pfile);

    for (int i=0; i<n*m; i++) {
        fwrite(&x[i],sizeof(double),1,pfile);
        fwrite(&y[i],sizeof(double),1,pfile);
    }
    fclose(pfile);
}
