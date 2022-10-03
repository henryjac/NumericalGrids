#include "grid_functions.h"

#include <memory>
#include <iostream>

std::array<std::shared_ptr<GFkt>,2> GFkt::pd() {
    std::array<std::shared_ptr<GFkt>,2> dxy;
    dxy[0] = std::make_shared<GFkt>(grid);
    dxy[1] = std::make_shared<GFkt>(grid);
    std::shared_ptr<Matrix> pdij;
    for (int i=0; i<grid->xsize(); i++) {
        for (int j=0; j<grid->ysize(); j++) {
            pdij = pd(i,j);
            dxy[0]->u(i,j) = (*pdij)(0,0);
            dxy[1]->u(i,j) = (*pdij)(1,0);
        }
    }
    return dxy;
}

GFkt GFkt::pdx() {
    GFkt dx(grid);
    for (int i=0; i<grid->xsize(); i++) {
        for (int j=0; j<grid->ysize(); j++) {
            dx.u(i,j) = (*pd(i,j))(0,0);
        }
    }
    return dx;
}

GFkt GFkt::pdy() {
    GFkt dy(grid);
    for (int i=0; i<grid->xsize(); i++) {
        for (int j=0; j<grid->ysize(); j++) {
            dy.u(i,j) = (*pd(i,j))(1,0);
        }
    }
    return dy;
}

// Two ways to calculate the laplacian, gives same result
GFkt GFkt::laplace() {
    GFkt laplace(grid);
    std::shared_ptr<Matrix> pdd;
    for (int i=0; i<grid->xsize(); i++) {
        for(int j=0; j<grid->ysize(); j++) {
            pdd = pd2(i,j);
            laplace.u(i,j) = (*pdd)(0,0) + (*pdd)(2,0);
        }
    }
    return laplace;
}

GFkt GFkt::laplace2() {
    GFkt x(grid);
    GFkt y(grid);
    x = pdx();
    y = pdy();
    x = x.pdx();
    y = y.pdy();
    return x + y;
}

DiffCase GFkt::getDiffCase(int i, int j) {
    if (i == 0 && j == 0) return BOTTOM_LEFT;
    if (i == 0 && j == grid->ysize()-1) return TOP_LEFT;
    if (i == grid->xsize()-1 && j == 0 ) return BOTTOM_RIGHT;
    if (i == grid->xsize()-1 && j == grid->ysize()-1) return TOP_RIGHT;
    if (i == 0) return LEFT;
    if (i == grid->xsize()-1) return RIGHT;
    if (j == 0) return BOTTOM;
    if (j == grid->ysize()-1) return TOP;
    return INSIDE;
}

void set_indices(bool change, int start_index, int inds[], int length) {
    if (change) {
        for (int i=0; i<length; i++) {
            // multiply by -1 for odd indices
            inds[i] = start_index + (2*((i+1)%2)-1)*(i/2+1);
        }
    } else {
        for (int i=0; i<length; i++) {
            inds[i] = start_index;
        }
    }
}

double one_sided_diff(int factor, double s, double s1, double s2, double h) {
    // factor should be 1 for left sided difference, -1 for right sided
    return factor*(3*s - 4*s1 + s2) / (2*h);
}

double central_diff(double sp1, double sm1, double h) {
    return (sp1 - sm1) / (2*h);
}

template <class intFnc>
double one_sided_diff_2d(intFnc f, int fc1, int fc2, int is[3], int js[3], double hs[2]) {
    double s[3];
    for (int i=0; i<3; i++) {
        s[i] = one_sided_diff(fc1, f(is[i],js[0]), f(is[i],js[1]), f(is[i],js[2]), hs[1]);
    }
    return one_sided_diff(fc2, s[0], s[1], s[2], hs[0]);
}

template <class intFnc>
double one_sided_central_diff(intFnc f, DiffDirection1D dir,  int fc, int is[3], int js[3], double hs[2]) {
    double s[3];
    double h;
    switch(dir) {
        case(XI): {
            h = hs[1];
            for (int i=0; i<3; i++)
                s[i] = central_diff(f(is[i],js[1]), f(is[i],js[2]), hs[0]);
            break;
        } case(ETA): {
            h = hs[0];
            for (int j=0; j<3; j++) 
                s[j] = central_diff(f(is[1],js[j]), f(is[2],js[j]), hs[1]);
            break;
        }
    }
    return one_sided_diff(fc, s[0], s[1], s[2], h);
}

// Implemented the same way as for second derivatives,
// not as efficient maybe, but generalizes the concept
std::shared_ptr<Matrix> GFkt::pd(int i, int j) {
    auto x = [&](int i, int j) -> double {return grid->getx(i,j);};
    auto y = [&](int i, int j) -> double {return grid->gety(i,j);};
    auto f = [&](int i, int j) -> double {return u.operator()(i,j);};

    // e: xi, n: eta (how the greek letters look)
    double xe = pderiv(x,XI,i,j), xn = pderiv(x,ETA,i,j);
    double ye = pderiv(y,XI,i,j), yn = pderiv(y,ETA,i,j);

    double ue = pderiv(f, XI, i, j);
    double un = pderiv(f, ETA, i, j);

    SquareMatrix A(2);
    double Avals[] = {
        xe, ye, xn, yn
    };
    A.fillMatrix(Avals);

    Matrix b(3,1);
    double bvals[] = {
        ue, un
    };
    b.fillMatrix(bvals);

    std::shared_ptr<Matrix> X = A.SolveEq(b);
    return X;
}

template <class intFnc>
double GFkt::pderiv(intFnc f, DiffDirection1D curr_dir, int i, int j) {
    // f is a function defined on the grid
    double h;
    int curr_ind, size;
    int ipm[4];
    int jpm[4];
    switch(curr_dir) {
        case(XI): {
            curr_ind = i; size = grid->xsize();
            set_indices(true,i,ipm,4);
            set_indices(false,j,jpm,4);
            h = h_xi;
            break;
        }
        case(ETA): {
            curr_ind = j; size = grid->ysize();
            set_indices(false,i,ipm,4);
            set_indices(true,j,jpm,4);
            h = h_eta;
            break;
        }
    }
    if (curr_ind == 0) {
        return one_sided_diff(-1, f(i,j), f(ipm[0],jpm[0]), f(ipm[2],jpm[2]), h);
    } else if (curr_ind == size-1) {
        return one_sided_diff(1, f(i,j), f(ipm[1],jpm[1]), f(ipm[3],jpm[3]), h);
    } else {
        return central_diff(f(ipm[0],jpm[0]), f(ipm[1],jpm[1]), h);
    }
}

std::shared_ptr<Matrix> GFkt::pd2(int i, int j) {
    auto x = [&](int i, int j) -> double {return grid->getx(i,j);};
    auto y = [&](int i, int j) -> double {return grid->gety(i,j);};
    auto f = [&](int i, int j) -> double {return u.operator()(i,j);};

    double xee = pderiv2(x,XI_XI,i,j), xnn = pderiv2(x,ETA_ETA,i,j);
    double yee = pderiv2(y,XI_XI,i,j), ynn = pderiv2(y,ETA_ETA,i,j);
    double xen = pderiv2(x,XI_ETA,i,j), yen = pderiv2(y,XI_ETA,i,j);

    double uee = pderiv2(f,XI_XI,i,j), unn = pderiv2(f,ETA_ETA,i,j);
    double uen = pderiv2(f,XI_ETA,i,j);

    double xe = pderiv(x,XI,i,j), xn = pderiv(x,ETA,i,j);
    double ye = pderiv(y,XI,i,j), yn = pderiv(y,ETA,i,j);

    std::shared_ptr<Matrix>uxy = pd(i,j);
    double ux = (*uxy)(0,0);
    double uy = (*uxy)(1,0);

    SquareMatrix A(3);
    double Avals[] = {
        xe*xe, 2*xe*ye,     ye*ye,
        xe*xn, xe*yn+xn*ye, ye*yn,
        xn*xn, 2*xn*yn,     yn*yn,
    };
    A.fillMatrix(Avals);

    Matrix b(3,1);
    double bvals[] = {
        uee - ux*xee - uy*yee,
        uen - ux*xen - uy*yen,
        unn - ux*xnn - uy*ynn,
    };
    b.fillMatrix(bvals);
    
    std::shared_ptr<Matrix> X = A.SolveEq(b);
    return X;
}

template <class intFnc>
double GFkt::pderiv2(intFnc f, DiffDirection2D curr_diff_dir, int i, int j) {
    DiffCase curr_diff_case = getDiffCase(i, j);
    double h;
    int curr_ind, size;
    int ipm[6];
    int jpm[6];
    switch (curr_diff_dir) {
        case(XI_XI) : { 
            curr_ind = i;
            size = grid->xsize();
            set_indices(true,i,ipm,6);
            set_indices(false,j,jpm,6);
            h = h_xi;
            break;
        } case(ETA_ETA) : {
            curr_ind = j;
            size = grid->ysize();
            set_indices(false,i,ipm,6);
            set_indices(true,j,jpm,6);
            h = h_eta;
            break;
        } case(XI_ETA) : {
            double hs[2] = {h_xi, h_eta};
            set_indices(true,i,ipm,6);
            set_indices(true,j,jpm,6);
            switch(curr_diff_case) {
                case BOTTOM_LEFT : {
                    int is[3] = {i,ipm[0],ipm[2]}, js[3] = {j,jpm[0],jpm[2]};
                    return one_sided_diff_2d(f,-1,-1,is,js,hs);
                }
                case TOP_LEFT : {
                    int is[3] = {i,ipm[0],ipm[2]}, js[3] = {j,jpm[1],jpm[3]};
                    return one_sided_diff_2d(f,1,-1,is,js,hs);
                }
                case BOTTOM_RIGHT : {
                    int is[3] = {i,ipm[1],ipm[3]}, js[3] = {j,jpm[0],jpm[2]};
                    return one_sided_diff_2d(f,-1,1,is,js,hs);
                }
                case TOP_RIGHT : {
                    int is[3] = {i,ipm[1],ipm[3]}, js[3] = {j,jpm[1],jpm[3]};
                    return one_sided_diff_2d(f,-1,-1,is,js,hs);
                }
                case LEFT : {
                    int is[3] = {i,ipm[0],ipm[2]}, js[3] = {j,jpm[0],jpm[1]};
                    return one_sided_central_diff(f,XI,-1,is,js,hs);
                }
                case RIGHT : {
                    int is[3] = {i,ipm[1],ipm[3]}, js[3] = {j,jpm[0],jpm[1]};
                    return one_sided_central_diff(f,XI,1,is,js,hs);
                }
                case BOTTOM : {
                    int is[3] = {i,ipm[0],ipm[1]}, js[3] = {j,jpm[0],jpm[2]};
                    return one_sided_central_diff(f,ETA,-1,is,js,hs);
                }
                case TOP : {
                    int is[3] = {i,ipm[0],ipm[1]}, js[3] = {j,jpm[1],jpm[3]};
                    return one_sided_central_diff(f,ETA,1,is,js,hs);
                }
                default: {
                    return (f(ipm[0],jpm[0]) - f(ipm[0],jpm[1]) - f(ipm[1],jpm[0]) + f(ipm[1],jpm[1])) / 
                        (4*h_xi*h_eta);
                }
            }
        }
    }
    double ret;
    if (curr_ind == 0) {
        ret = (2*f(i,j) - 5*f(ipm[0],jpm[0]) + 4*f(ipm[2],jpm[2]) - f(ipm[4],jpm[4])) / (h*h);
    } else if (curr_ind == size-1) {
        ret = (2*f(i,j) - 5*f(ipm[1],jpm[1]) + 4*f(ipm[3],jpm[3]) - f(ipm[5],jpm[5])) / (h*h);
    } else {
        ret = (f(ipm[0],jpm[0]) - 2*f(i,j) + f(ipm[1],jpm[1])) / (h*h);
    }
    return ret;
}
