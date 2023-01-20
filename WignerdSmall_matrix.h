#include <iostream>
#include "Eigen/Core"
#define PRINT_MAT(X) cout << #X << ":\n" << X << endl << endl   //行列の出力

Eigen::MatrixXd wigner_d_small_pihalf(int dimension) {
    Eigen::MatrixXd d;
    switch (dimension) {
        case 0:
            d << 1;
            return d;
            break;
        
        case 1:
            d <<       0.5,  0.707107,       0.5,
                 -0.707107,         0,  0.707107,
                       0.5, -0.707107,       0.5;
            return d;
            break;
    }
 }

 int main () {
    wigner_d_small_pihalf
 }