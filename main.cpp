#include <iostream>
#include <vector>
#include <cstdlib>  // abs() for integer
#include <numeric>  //accumulate
#include <cmath>
#include <time.h>
#include <omp.h> 
#include "MatrixElement.h"  //Hamiltonian matrix element
#include "Eigenstate.h" //Calculation of eigen vector components
#include "diagonalization_dsyevd.h" //DSYEVD diagonalizing calculation

#define PRINT_MAT(X) cout << #X << ":\n" << X << endl << endl   //行列の出力
#define EIGEN_DONT_PARALLELIZE
using namespace std;

int main() 
{   
    //time log
    ofstream time("result_run_time.txt"); 
    time << "n0 = " << n0 << "\t" << "npls = " << npls << "\t" << "nmns = " << nmns << "\t" << "j = " << j << "\t" << "m = " << calcm << "\n" ;
    time.close();
    clock_t start_whole = clock();
    
    int m = calcm;

    MatrixParameter hoge(m, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
    int dimension = hoge.matrix_size();
    cout << dimension << endl;

    //Hamiltonian
    clock_t start_matrix = clock();
    //Eigen::MatrixXd Hamiltonian =  Ts + Vs + Tb + Vb + Vbb + Tr + Vr + Vsb + Vsr + Vbr + Vam + VamR;
    Eigen::MatrixXd Hamiltonian = H_H(m);
    clock_t end_matrix = clock();
    const double time_matrix = static_cast<double>(end_matrix - start_matrix) / CLOCKS_PER_SEC;
    printf("matrix time %lf[s]\n", time_matrix);
 
    //diagonalizing
    clock_t start_diag = clock();
    Eigen::VectorXd eig_val = Eigen::VectorXd::Zero(Hamiltonian.rows());
    Eigen::MatrixXd eig_vec = Eigen::MatrixXd::Zero(Hamiltonian.rows(), Hamiltonian.rows());
    tie(eig_val, eig_vec) = diagonalization_dsyevd(Hamiltonian);    //diagonalization, 基底対称化を行わない場合, 実対称行列 from diagonalization_dsyevd.h
    eigen_state(m, eig_val, eig_vec);     //固有べベクトルの成分計算 from Eigenstate.h
    
    clock_t end_diag = clock();
    const double time_diag = static_cast<double>(end_diag - start_diag) / CLOCKS_PER_SEC;
    printf("diagonalization time %lf[s]\n", time_diag);
    
    clock_t end_whole = clock();
    const double time_whole = static_cast<double>(end_whole - start_whole) / CLOCKS_PER_SEC;
    printf("whole time %lf[s]\n", time_whole);

    //time log
    ofstream time2("result_run_time.txt", ios::app);
    time2 << "matrix time " << time_matrix << "[s]\n";
    time2 << "diagonalization time " << time_diag << "[s]\n";
    time2 << "whole time " << time_whole << "[s]\n";
    time2.close();

    return 0;
}