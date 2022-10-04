#include <iostream>
#include <exception>
#include <string.h>
#include <memory>

class Matrix {
    protected:
        int n,m; // rows / columns
        double *X;
    public:
        // Constructors & Destructor
        Matrix();
        Matrix(int,int);
        Matrix(const Matrix&);
        Matrix(Matrix&&);
        ~Matrix();
        // Overloading operators
        Matrix& operator=(const Matrix&);
        Matrix& operator+=(const Matrix&);
        Matrix& operator*=(const Matrix&);
        Matrix& operator*=(const double);
        Matrix& operator/=(const double);
        const Matrix operator*(const Matrix&) const;
        const Matrix operator+(const Matrix&) const;
        const Matrix operator-(const Matrix&) const; 
        const Matrix operator*(const double) const;
        const Matrix operator/(const double) const;
        double& operator()(int,int) const;
        // Other methods
        double norm_1() const;
        double norm_inf() const;
        void printMatrix() const;
        void fillMatrix(double*);
        double* get_vals() const;
        void set_val(double,int,int);
        int* size() const;
        int cols() const;
        int rows() const;
};

class SquareMatrix : public Matrix {
    public:
        SquareMatrix(int); // Empty matrix of size int
        SquareMatrix(const Matrix&); // Matrix copy-constructed
        SquareMatrix(int,int); // Diagonal matrix with diag-value
        SquareMatrix& pow(int);
        SquareMatrix& exp(double tol=1e-10);
        std::shared_ptr<Matrix> SolveEq(const Matrix&);
    private:
        // For solving linear systems, used partially
        // in project 4 (for 3x3 only), but I found it fun to implement it generally
        // Methods defined in project4
        std::array<std::shared_ptr<SquareMatrix>,2> LUdecomp();
        std::shared_ptr<Matrix> forward_sub(const Matrix&);
        std::shared_ptr<Matrix> backward_sub(const Matrix&);
};

class SizeError : public std::exception {
    private:
        int n1,m1,n2,m2;
        std::string operation;
        virtual const char* what() const throw() {
            std::string wrong;
            wrong.append("Invalid operation ");
            wrong.append(operation);
            wrong.append(" for matrices of size ");
            wrong.append(std::to_string(n1));
            wrong.append("x");
            wrong.append(std::to_string(m1));
            wrong.append(" and ");
            wrong.append(std::to_string(n2));
            wrong.append("x");
            wrong.append(std::to_string(m2));
            static char* arr = new char[wrong.length() + 1];
            strcpy(arr, wrong.c_str());
            return arr;
        }
    public :
        SizeError(const Matrix& A, const Matrix& B, std::string oper="+") {
            int *sizeA = A.size(); int *sizeB = B.size();
            n1 = *sizeA; m1 = *(sizeA+1);
            n2 = *sizeB; m2 = *(sizeB+1);
            operation = oper;
        }
};
