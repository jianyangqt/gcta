#include "Matrix.hpp"
#include "Eigen/Eigen"
#include <iostream>

using namespace std;
using namespace Eigen;

int main(int argc, char* argv[]){
    int n = 100;
    MatrixXd A(n, n);
    A.setRandom();
    A = A * A.adjoint();
    INVmethod method = INV_QR;
    double logdet = 0;
    int rank = 0;
    SquareMatrixInverse(A, logdet, rank, method);
    cout << method << endl;
    cout << logdet << endl;

    MatrixXd B(2,2);
    B << 1,1,
      1,1;
    INVmethod method2 = INV_FQR;
    SquareMatrixInverse(B, logdet, rank, method2);
    cout << "logdet:" << logdet << std::endl;
    cout << "rank: " << rank  << std::endl;
    cout << B << endl;
}
