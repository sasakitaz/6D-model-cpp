/*****************************************************************************************************
j_y行列の対角化を用いたアルゴリズム
参考: X. M. Feng, P. Wang, W. Yang, and G. R. Jin, Phys. Rev. E 92, 043307 – Published 19 October 2015
https://journals.aps.org/pre/abstract/10.1103/PhysRevE.92.043307

ちなみにアルゴリズムは違うがpythonでのライブラリ: 素直に階乗を計算するアルゴリズム
https://docs.sympy.org/latest/modules/physics/wigner.html

量子数j, mは整数の場合のみ実装
半整数の場合は未実装
*****************************************************************************************************/
#pragma once
#include <iostream>
#include <complex>
#include <cstdlib>  // abs() for integer
#include <math.h>
#include "tools.h"  //Wigner3jSymbol product, DVR, classの継承
#define PRINT_MAT(X) cout << #X << ":\n" << X << endl << endl   //行列の出力

using namespace std;

class WignerdSmall : public tools
{
protected:
    int p;   //spherical harmonics parameter
    int q;   //spherical harmonics parameter

public:
    //基底クラスのコンストラクタ呼び出し
    WignerdSmall(int temm, int temjrow, int temmBrow, int temkrow, int temn0row, int temvrow, int temlrow, int temjcolumn, int temmBcolumn, int temkcolumn, int temn0column, int temvcolumn, int temlcolumn) : tools(temm, temjrow, temmBrow, temkrow, temn0row, temvrow, temlrow, temjcolumn, temmBcolumn, temkcolumn, temn0column, temvcolumn, temlcolumn)
    {
    }
    double jy_element(int j_jy, int m_jy);  //j_zの固有状態基底によるj_y行列の行列要素
    Eigen::MatrixX<complex<double>> jy_eigen ();   //j_y行列の対角化
    Eigen::MatrixXd wigner_d_small(double beta);    //j_y行列の固有ベクトル(ユニタリー行列<j, m|j, mu>)を用いたd関数の複素Fourier展開の表式
};

double WignerdSmall::jy_element(int j_jy, int m_jy) {
    double Xm = sqrt((j_jy + m_jy)*(j_jy - m_jy +1));
    return Xm;
}

Eigen::MatrixX<complex<double>> WignerdSmall::jy_eigen () {
    complex<double> el1;
    complex<double> el2;
    //Eigen::MatrixX<complex<double>> jy(2*jrow + 1, 2*jrow + 1);
    Eigen::MatrixX<complex<double>> jy = Eigen::MatrixX<complex<double>>::Zero(2*jrow + 1, 2*jrow + 1);
    for (int i = 1; i <= 2*jrow; i++) {
        el1 = complex<double> (0, + jy_element(jrow, jrow - i + 1)/(-2)); //下三角
        el2 = complex<double> (0, - jy_element(jrow, - jrow + i)/(-2));   //上三角
        jy(i, i - 1) = el1;    //下三角
        jy(i - 1, i) = el2;    //上三角
    }
    Eigen::SelfAdjointEigenSolver< Eigen::MatrixX<complex<double>>> s(jy);   //三重対角行列なのでLAPACKにすれば高速化できる
    cout << "!!diagonalize_jy!!" << endl;
    return s.eigenvectors();
}

Eigen::MatrixXd WignerdSmall::wigner_d_small(double beta) {
    Eigen::MatrixX<complex<double>> jy_eigen_vec = jy_eigen();
    Eigen::MatrixXd d = Eigen::MatrixXd::Zero(2*jrow + 1, 2*jrow + 1);
    double dmn;
    complex<double> phase;
    complex<double> coeff;
    for (int m = - jrow; m <= jrow; m++) {
        for (int n = -jrow; n <= jrow; n++) {
            complex<double> dmn_comp = 0;
            for (int mu = - jrow; mu <= jrow; mu++){
                phase = complex<double> (cos(mu*beta), - sin(mu*beta));
                coeff = complex<double> (jy_eigen_vec(jrow + m, jrow + mu)*conj(jy_eigen_vec(jrow + n, jrow + mu)));
                dmn_comp += phase*coeff;
            //debag
            }
            if (abs(imag(dmn_comp)) > pow(10, -10)) {
                cout << "\t\t\t!!error!! dmn is complex" <<endl;
            }
            dmn = real(dmn_comp);
            if (abs(dmn) > pow(10, -10)) {
                d(jrow + m, jrow + n) = dmn;
            }
            else {
                d(jrow + m, jrow + n) = 0;
            }
        }
    }
    return d;
}