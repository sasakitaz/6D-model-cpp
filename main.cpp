#include <iostream>
#include <vector>
#include <cstdlib>  // abs() for integer
#include <numeric>  //accumulate
#include <cmath>
#include <time.h>
#include "MatrixElement.h"  //Hamiltonian項の行列要素
#include "Eigenstate.h" //固有値計算後の成分計算
#include "diagonalization_zheevr.h"
#include "diagonalization_dsyevd.h"
#include "SymmetrizeDiagonalization.h"  //基底対称化した行列の対角化

#define PRINT_MAT(X) cout << #X << ":\n" << X << endl << endl   //行列の出力

using namespace std;

int main() 
{   
    //time log
    ofstream time("result_run_time.txt");    //txtファイル書き出し
    time << "n0 = " << n0 << "\t" << "npls = " << npls << "\t" << "nmns = " << nmns << "\t" << "j = " << j << "\t" << "m = " << calcm << "\n" ;
    time.close();
    clock_t start_whole = clock();
    
    //Hamiltonian項の行列をMatrixElement.hから取り出す
    int m = calcm;

    //対称化0, 1, 3, 6の選択を標準入力にする
    int sym_size;
    cout << "input the method of basis symmetrization, 0, 1, 3, 6"<<endl;
    cout << "0. no symmetrize, dsyevd, default"<<endl;
    cout << "1. symmetrize and diagonarize unitary block matrix, for test"<<endl;
    cout << "3. symmetrize and diagnarize A, E, F matrix，if eigenvectors are required, zheevr (under construction)"<<endl;
    cout << "6. symmetrize and diagnarize A1, E1, F1 matrix，if eigenvectors are not required, zheevr (under construction)"<<endl;
    cin >> sym_size;

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

    if (sym_size == 0) {diagonalization (m, Hamiltonian);}
    else if (sym_size == 1) {symmetrized_diagonalization1 (m, Hamiltonian);}
    else if (sym_size == 3) {symmetrized_diagonalization3 (m, Hamiltonian);}
    else if (sym_size == 6) {symmetrized_diagonalization6 (m, Hamiltonian);}
    else {cout << "ERROR sym_size out of range"<<endl;}


    clock_t end_diag = clock();
    const double time_diag = static_cast<double>(end_diag - start_diag) / CLOCKS_PER_SEC;
    printf("diagonalization time %lf[s]\n", time_diag);
    
    clock_t end_whole = clock();
    const double time_whole = static_cast<double>(end_whole - start_whole) / CLOCKS_PER_SEC;
    printf("whole time %lf[s]\n", time_whole);

    //time log
    ofstream time2("result_run_time.txt", ios::app);    //txtファイル書き出し
    time2 << "matrix time " << time_matrix << "[s]\n";
    time2 << "diagonalization time " << time_diag << "[s]\n";
    time2 << "whole time " << time_whole << "[s]\n";
    time2.close();

    return 0;
}