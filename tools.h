#pragma once
#include <iostream>
#include <vector>
#include <cstdlib>  // abs() for integer
#include <numeric>  //accumulate
#include <cmath>
#include <bits/stdc++.h>    //pi: M_PI
#include "Wigner3jSymbol.h" //wigner 3j symbol
#include "Parameter.h"
#include "Eigen/Core"  
#include "Eigen/Dense"  
#define PRINT_MAT(X) cout << #X << ":\n" << X << endl << endl 
using namespace std;


class tools : public MatrixParameter
{
protected:
    double kzz;     //Morse potential parameter
    double az;      //Morse potential parameter 
    int p;          //spherical harmonics parameter
    int q;          //spherical harmonics parameter
    int m_couple;   //angular momentum coupling in Vam and VamR
public:
    tools(int temm, int temjrow, int temmBrow, int temkrow, int temn0row, int temvrow, int temjcolumn, int temmBcolumn, int temkcolumn, int temn0column, int temvcolumn) : MatrixParameter(temm, temjrow, temmBrow, temkrow, temn0row, temvrow, temjcolumn, temmBcolumn, temkcolumn, temn0column, temvcolumn)
    {
    }
    double Wigner3jSymbol(int p, int q);
    double Wigner3jSymbol_couple(int p, int q, int m_couple);
};

double tools::Wigner3jSymbol(int p, int q) {
    double symbol1 = wigner_3j(    jrow, p,  jcolumn, 
                                 - krow, q,  kcolumn);
    double symbol2 = wigner_3j(    jrow, p,  jcolumn,
                                - mBrow, 0, mBcolumn);
    double symbol = symbol1*symbol2;
    return symbol;
}

double tools::Wigner3jSymbol_couple(int p, int q, int m_couple) {
    double symbolpls1 = wigner_3j(    jrow,        p,  jcolumn, 
                                    - krow,        q,  kcolumn);
    double symbolpls2 = wigner_3j(    jrow,        p,  jcolumn,
                                   - mBrow, m_couple, mBcolumn);
    double symbolpls = symbolpls1*symbolpls2;
    return symbolpls;
}

tuple<Eigen::MatrixXd, Eigen::MatrixXd>  DVR(double kzz, double az) {
    Eigen::MatrixXd x = Eigen::MatrixXd::Zero(n0 + 1, n0 + 1);
    for (int i = 1; i <= n0; i++) {
        x(i, i - 1) = (double)sqrt(i);    //creation operator
        x(i - 1, i) = (double)sqrt(i);    //annihilation operator
    }
    x = sqrt(shbaromega/(4*kzz*az*az))*x;   //position operator
    Eigen::SelfAdjointEigenSolver< Eigen::Matrix<double, n0 + 1, n0 + 1>> s(x); 
    cout << "!!diagonalize_DVR!!" << endl;
    return forward_as_tuple(s.eigenvalues(), s.eigenvectors());
}  

Eigen::MatrixXd getRi() {
    Eigen::MatrixXd R;
    tie(R, ignore) = DVR(kzz, az);
    return R;
}

Eigen::MatrixXd getRi_vec() {
    Eigen::MatrixXd R_vec;
    tie(ignore, R_vec) = DVR(kzz, az);
    return R_vec;
}

Eigen::MatrixXd Ri = getRi();
Eigen::MatrixXd Ri_vec = getRi_vec();

