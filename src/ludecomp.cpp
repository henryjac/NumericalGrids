#include "matrix.h"

#include <memory>
#include <array>

std::shared_ptr<Matrix> SquareMatrix::SolveEq(const Matrix& b) {
    std::array<std::shared_ptr<SquareMatrix>,2> LU = LUdecomp();
    std::shared_ptr<Matrix> y = LU[0]->forward_sub(b);
    std::shared_ptr<Matrix> x = LU[1]->backward_sub(*y);
    return x;
}

// TODO: Implement partial pivoting
std::array<std::shared_ptr<SquareMatrix>,2> SquareMatrix::LUdecomp() {
    std::array<std::shared_ptr<SquareMatrix>, 2> LU;
    for (int i=0; i<2; i++)
        LU[i] = std::make_shared<SquareMatrix>(n);

    std::shared_ptr<SquareMatrix> P = std::make_shared<SquareMatrix>(n,1);

    for (int i=0; i<n; i++) {
        for (int j=i; j<n; j++) {
            double sum(0);
            for (int k=0; k<i; k++) {
                sum += LU[0]->X[i*n+k] * LU[1]->X[k*n+j];
            }
            LU[1]->X[i*n+j] = X[i*n+j] - sum;
        }
        for (int j=i; j<n; j++) {
            if (i==j) {
                LU[0]->X[i*n+i] = 1;
            } else {
                double sum(0);
                for (int k=0; k<i; k++) {
                    sum += LU[0]->X[j*n+k] * LU[1]->X[k*n+i];
                }
                LU[0]->X[j*n+i] = (X[j*n+i] - sum) / LU[1]->X[i*n+i];
            }
        }
    }
    return LU;
}

std::shared_ptr<Matrix> SquareMatrix::forward_sub(const Matrix& b) {
    std::shared_ptr<Matrix> x = std::make_shared<Matrix>(n,1);
    for (int i=0; i<n; i++) {
        double sum(0);
        for (int j=0; j<i; j++) {
            sum += X[i*n+j] * (*x)(j,0);
        }
        (*x)(i,0) = (b(i,0) - sum) / X[i*n+i];
    }
    return x;
}

std::shared_ptr<Matrix> SquareMatrix::backward_sub(const Matrix& b) {
    std::shared_ptr<Matrix> x = std::make_shared<Matrix>(n,1);
    for (int i=n-1; i>=0; i--) {
        double sum(0);
        for (int j=i+1; j<n; j++) {
            sum += X[i*n+j] * (*x)(j,0);
        }
        (*x)(i,0) = (b(i,0) - sum) / X[i*n+i];
    }
    return x;
}
