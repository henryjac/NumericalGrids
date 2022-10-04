#include "matrix.h"

#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <exception>

// Matrix: method definitions
Matrix::Matrix() : n(), m(), X(nullptr) {};

// (n,m)-matrix with zero entries
Matrix::Matrix(int nn, int mm) : n(nn), m(mm) {
    X = new double[n*m];
    std::fill(X,X+m*n,0.0);
}
// Copy constructor
Matrix::Matrix(const Matrix& P) : n(P.n), m(P.m) {
    X = new double[n*m];
    std::copy(P.X,P.X+n*m,X);
}
// Move assignement
Matrix::Matrix(Matrix&& P) : n(P.n), m(P.m), X(P.X) {
    P.m = 0; P.n = 0; P.X = nullptr;
}
Matrix::~Matrix() {delete[] X;} // Delete dynamically created array

// Copy operator
Matrix& Matrix::operator=(const Matrix& P) {
    if (this != &P) {
        if (n*m != P.n*P.m) {
            if (X != nullptr) {
                delete[] X;
            }
            // Create from P if only it isn't 0x0-matrix
            if (P.X != nullptr) {
                X = new double[P.n*P.m];
            }
        }
        n = P.n; 
        m = P.m; 
        // Use std lib for copying the array
        std::copy(P.X,P.X+n*m,X);
    }
    return *this; 
};

Matrix& Matrix::operator+=(const Matrix& P) {
    // Check if matrix sizes are compatible for addition
    try {
        if (n != P.n && m != P.m) {SizeError size_err(*this, P, "+"); throw size_err;}
    } catch (std::exception& e) {
        std::cerr << e.what() << std::endl; 
        exit(1);
    }
    for (int i=0; i<n*m; i++) {
        this->X[i] += P.X[i];
    }
    return *this;
}

Matrix& Matrix::operator*=(const Matrix& P) {
    // Check if matrix sizes are compatible for multiplication
    try {
        if (m != P.n) {SizeError size_err(*this, P, "*"); throw size_err;}
    } catch (std::exception& e) {
        std::cerr << e.what() << std::endl;
        exit(1);
    }
    Matrix newMatrix(n,P.m);
    for (int i=0; i<n; i++) {
        for (int j=0; j<P.m; j++) {
            for (int k=0; k<m; k++) {
                newMatrix.X[i*m+j] += X[i*m+k] * P.X[j+k*P.m];
            }
        }
    }
    *this = newMatrix;
    return *this;
}

Matrix& Matrix::operator*=(const double a) {
    for (int i=0; i<n*m; i++) {
        X[i]*=a;
    }
    return *this;
}

Matrix& Matrix::operator/=(const double a) {
    for (int i=0; i<n*m; i++) {
        X[i]/=a;
    }
    return *this;
}

const Matrix Matrix::operator*(const Matrix& P) const {
    // Check if matrix sizes are compatible for multiplication
    try {
        if (m != P.n) {SizeError size_err(*this, P, "*"); throw size_err;}
    } catch (std::exception& e) {
        std::cerr << e.what() << std::endl;
        exit(1);
    }
    Matrix newMatrix(n,P.m);
    for (int i=0; i<n; i++) {
        for (int j=0; j<P.m; j++) {
            for (int k=0; k<m; k++) {
                newMatrix.X[i+j*n] += X[i*n+k] * P.X[k+j*P.n];
            }
        }
    }
    return newMatrix;
}

const Matrix Matrix::operator+(const Matrix& P) const {
    Matrix newMatrix(n,m);
    for (int i=0; i<n*m; i++) {
        newMatrix.X[i] = this->X[i] + P.X[i];
    }
    return newMatrix;
}

const Matrix Matrix::operator-(const Matrix& P) const {
    Matrix newMatrix(n,m);
    for (int i=0; i<n*m; i++) {
        newMatrix.X[i] = this->X[i] - P.X[i];
    }
    return newMatrix;
}


const Matrix Matrix::operator*(const double a) const {
    Matrix newMatrix(n,m);
    for (int i=0; i<n*m; i++) {
        newMatrix.X[i] = this->X[i] * a;
    }
    return newMatrix;
}

// Doesnt seem to work as expected
const Matrix Matrix::operator/(const double a) const {
    Matrix newMatrix(n,m);
    for (int i=0; i<n*m; i++) {
        newMatrix.X[i] = this->X[i]/a;
    }
    return newMatrix;
}

double& Matrix::operator()(int i, int j) const {return X[i*m+j];}

// Calculates the 1-norm of the matrix
double Matrix::norm_1() const {
    double max(0);
    for (int j=0; j<m; j++) {
        double colSum(0);
        for (int i=0; i<n; i++) {
            colSum += abs(X[i*m+j]);
        }
        if (max < colSum) {
            max = colSum;
        }
    }
    return max;
}

// Calculates the inf-norm of the matrix
double Matrix::norm_inf() const {
    double max(0);
    for (int i=0; i<n; i++) {
        double rowSum(0);
        for (int j=0; j<m; j++) {
            rowSum += abs(X[i*m+j]);
        }
        if (max < rowSum) {
            max = rowSum;
        }
    }
    return max;
}

// Prints the matrix in classic format with hard edges around the elements
void Matrix::printMatrix() const {
    if (X == nullptr) {return;}
    std::cout << " _";
    for (int j=0; j<m*7+m; j++) {
        std::cout << " ";
    }
    std::cout << " _\n";
    int len=0;
    for (int i=0; i<n; i++) {
        std::cout << "\uFF5C ";
        std::cout << std::fixed;
        for (int j=0; j<m; j++) {
            for (float countdown=abs(X[i*m+j]); countdown>=10; countdown/=10) {
                len++;
            }
            // If number is bigger than 1e6 use scientific notation
            if (len > 6) {
                std::cout << std::scientific;
                /* std::cout << std::setprecision(1); */
                if (X[i*m+j] < 0) {
                    std::cout << std::setprecision(0);
                } else {
                    std::cout << std::setprecision(1);
                }
            }
            else {
                std::cout << std::fixed;
                if (X[i*m+j] < 0) {
                    std::cout << std::setprecision(4-len);
                } else {
                    std::cout << std::setprecision(5-len);
                }
            }
            std::cout << X[i*m+j] << " ";
            len=0;                                             
            std::cout << std::fixed;
        }
        std::cout << "\uFF5c\n";
    }
    std::cout << " \u203E";
    for (int i=0; i<m*7+m+1; i++) std::cout << " ";
    std::cout << "\u203E\n";
    std::cout << std::scientific;
}

void Matrix::fillMatrix(double vals[]) {
    // Copy array to X
    int size = n*m;
    for (int i=0; i<size; i++) {
        *(X+i) = *(vals+i);
    }
};

double* Matrix::get_vals() const {
    return X;
}

void Matrix::set_val(double x, int i, int j) {
    X[i*n+j] = x;
}

int* Matrix::size() const {
    int *size = new int[2];
    size[0] = n; size[1] = m;
    return size;
}

int Matrix::rows() const {
    return n;
}
int Matrix::cols() const {
    return m;
}

// SquareMatrix: method definitions
SquareMatrix::SquareMatrix(int nn) : Matrix(nn,nn) {};
SquareMatrix::SquareMatrix(const Matrix& P) : Matrix(P.rows(),P.cols()) {
    n = P.rows(); m = P.cols();
    if (n != m) {
        std::cout << "Can't make " << n << "x" << n << " matrix square." << std::endl;
        exit(1);
    }
    X = P.get_vals();
};

SquareMatrix::SquareMatrix(int nn, int x) : Matrix(nn,nn) {
    for (int i=0; i<n; i++) {
        X[i+i*n] = x;
    }
};

SquareMatrix& SquareMatrix::pow(int p) {
    if (p == 0) {
        SquareMatrix *idMatrix = new SquareMatrix(n,1);
        return *idMatrix;
    } else if (p < 0) {
        std::cout << "Matrix inversion not implemented." << std::endl;
        exit(1);
    }
    SquareMatrix *powMatrix = new SquareMatrix(*this);
    for (int i=0; i<p-1; i++) {
        *powMatrix *= *this;
    }
    return *powMatrix;
}

SquareMatrix& SquareMatrix::exp(double tol) {
    SquareMatrix *expMatrix = new SquareMatrix(n,1);
    SquareMatrix powMatrix(n,1);
    double err;
    int n=1;
    while (1) {
        powMatrix *= *this;
        powMatrix /= n;
        err = powMatrix.norm_1();
        if (err < tol) break;
        *expMatrix += powMatrix;
        n++;
    }
    return *expMatrix;
}
