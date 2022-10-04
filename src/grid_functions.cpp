#include "grid_functions.h"

#include <memory>
#include <iostream>
#include <cstdio>

GFkt::GFkt(std::shared_ptr<Domain> grid_) : u(grid_->xsize(),grid_->ysize()), grid(grid_) {
    if (grid->xsize() == 0 || grid->ysize() == 0) {
        std::cerr << "Grid must be instantiated before construction of GFkt object." << std::endl;
        exit(1);
    }
    h_xi = (double)1/(grid->xsize()-1);
    h_eta = (double)1/(grid->ysize()-1);
}

GFkt& GFkt::operator=(const GFkt& U) {
    if (this != &U) {
        u = U.u;
        grid = U.grid; // Share grid
    }
    return *this;
}

GFkt GFkt::operator+(const GFkt& U) const {
    if (grid == U.grid) {
        GFkt tmp(grid);
        tmp.u = u+U.u;
        return tmp;
    } else {
        std::cerr << "Grid functions defined on different grids." << std::endl;
        exit(1);
    }
}

GFkt GFkt::operator-(const GFkt& U) const {
    if (grid == U.grid) {
        GFkt tmp(grid);
        tmp.u = u-U.u;
        return tmp;
    } else {
        std::cerr << "Grid functions defined on different grids." << std::endl;
        exit(1);
    }
}

GFkt GFkt::operator*(const GFkt& U) const {
    if (grid == U.grid) {
        GFkt tmp(grid);
        for (int i=0; i<grid->xsize(); i++) {
            for (int j=0; j<grid->ysize(); j++) {
                tmp.u(i,j) = u(i,j) * U.u(i,j);
            }
        }
        return tmp;
    } else {
        std::cerr << "Grid functions defined on different grids." << std::endl;
        exit(1);
    }
}

GFkt GFkt::operator*(double a) const {
    GFkt tmp(grid);
    tmp.u = u*a;
    return tmp;
}

GFkt GFkt::operator/(double a) const {
    GFkt tmp(grid);
    tmp.u = u/a;
    return tmp;
}

void GFkt::fill_matrix(Fnc2D f) {
    for (int i=0; i<grid->xsize(); i++) {
        for (int j=0; j<grid->ysize(); j++) {
            u(i,j) = f(grid->operator()(i,j));
        }
    }
}

void GFkt::save(const char* valfile, const char* gridfile, const char* boundaryfile) {
    if (gridfile) {
        grid->save_grid(gridfile);
    }
    if (boundaryfile) {
        grid->save_boundary(boundaryfile);
    }
    FILE* pfile = fopen(valfile, "wb");
    // Write the shape so Julia can determine it explicitly
    int n = grid->xsize(), m = grid->ysize();
    fwrite(&n, sizeof(int), 1, pfile);
    fwrite(&m, sizeof(int), 1, pfile);

    double uij;
    for (int i=0; i<n; i++) {
        for (int j=0; j<m; j++) {
            uij = u(i,j);
            fwrite(&uij,sizeof(double),1,pfile);
        }
    }
    fclose(pfile);
}
