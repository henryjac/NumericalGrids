#include "../../project2/src/matrix.h"

#include <array>
#include <memory>
#include <iostream>

int main() {
    std::shared_ptr<SquareMatrix> A = std::make_shared<SquareMatrix>(3);
    double Avals[9] = {1,1,0,1,0,1,1,2,0};
    A->fillMatrix(Avals);

    Matrix b(3,1);
    double bvals[3] = {1,0,0};
    b.fillMatrix(bvals);

    std::shared_ptr<Matrix> x = A->SolveEq(b);
    std::cout << "A =\n";
    A->printMatrix();
    std::cout << "x =\n";
    x->printMatrix();
    std::cout << "b =\n";
    b.printMatrix();
    std::cout << "A\\b =\n";
    (*A**x).printMatrix();
}
