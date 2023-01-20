// zheev.cc
// 複素エルミート行列の固有値方程式を解く
#pragma once
#include <stdio.h>
#include <complex>
#include "MatrixElement.h"   // x(行, 列)
typedef std::complex<double> Complex;
#define PRINT_MAT(X) cout << #X << ":\n" << X << endl << endl   //行列の出力

extern "C" {
  void zheev_ ( const char& JOBZ, const char& UPLO,
		const int& N, Complex** A, const int& LDA,
		double* W, Complex* WORK, const int& LWORK,
		double* RWORK,
		int& INFO, int JOBZlen, int UPLOlen );
};

// 複素エルミート行列 H の固有値方程式 H x = E x を解く簡易関数
// Input: H[N][N],  Output: U[N][N], E[N]
// ただし、H[i][j<=i] の下三角しか参照されない
// U^\dagger H U が対角行列となる
//
template <int N> int zheev( Complex H[N][N], Complex U[N][N], double E[N] )
{
  int i, j, info;
  static Complex cwork[6*N];
  static double  rwork[3*N - 2];

  for( i=0; i<N; i++ ){
    for( j=0; j<N; j++ ){
      U[i][j] = H[j][i];
    }
  }

  zheev_( 'V', 'L', N, (Complex**)U, N, E, cwork, 6*N, rwork, info, 1, 1 );

  for( i=0; i<N; i++ ){
    for( j=i+1; j<N; j++ ){
      Complex temp;
      temp = U[i][j];
      U[i][j] = U[j][i];
      U[j][i] = temp;
    }
  }

  return info;
}

int diagonalization(Eigen::MatrixXd mat) {
    const int N = 9;
    Complex H[N][N], U[N][N];
    double E[N];
    int i, j, info;

    for (i = 0; i < N; i++) {
        for(j = 0; j <= i; j++) {
            H[i][j] = mat(i, j);
    }
  }
  
    printf("# H\n");
    for( i=0; i<N; i++ ){
        for( j=0; j<N; j++ ){
        printf("%+f ", real(H[i][j]) );
    }
    printf("\n");
    }

    info = zheev( H, U, E );

    printf("# info=%d.\n", info );

    printf("# Eigen values.\n");
    for( j=0; j<N; j++ ){
        printf("%+f ", E[j] );
    }
    printf("\n");

    printf("# Eigen vectors.\n");
    for( i=0; i<N; i++ ){
        for( j=0; j<N; j++ ){
        printf("%+f ", real(U[i][j]) );
    }
    printf("\n");
    }

    return 0;
}
