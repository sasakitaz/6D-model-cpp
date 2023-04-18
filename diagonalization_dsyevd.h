/***********************************************************************************************
2023.04.14 追記
free( (void*)a );
free( (void*)w );

2023.04.18 追記
ポインタを使って配列aへの移し替えを省略：a = mat.data();
************************************************************************************************/
#pragma once
#include <stdlib.h>
#include <stdio.h>
#include "Parameter.h"
/* DSYEVD prototype */
extern "C" {
    void dsyevd_( const char*  jobz, const char*  uplo, int* n, double* a, int* lda,
                    double* w, double* work, int* lwork, int* iwork, int* liwork, int* info );
}

/* Auxiliary routine: printing a matrix */
void print_matrix( const char*  desc, int m, int n, double* a, int lda ) {
        int i, j;
        printf( "\n %s\n", desc );
        for( i = 0; i < m; i++ ) {
                for( j = 0; j < n; j++ ) printf( " %6.2f", a[i+j*lda] );
                printf( "\n" );
        }
}

/* Auxiliary routines prototypes 
extern void print_matrix( const char&  desc, int m, int n, double* a, int lda );
*/

/* Main program */
tuple<Eigen::VectorXd, Eigen::MatrixXd> diagonalization_dsyevd(Eigen::MatrixXd mat) {
        /* Parameters */
        #define N mat.rows()
        #define LDA N
        /* Locals */
        int n = N, lda = LDA, info, lwork, liwork;
        int iwkopt;
        int* iwork;
        double wkopt;
        double* work;
        /* Local arrays */
        double w[N];
        double* a;
        a = mat.data();

        /* Executable statements */
        /* Query and allocate the optimal workspace */
        lwork = -1;
        liwork = -1;
        dsyevd_( "V", "L", &n, a, &lda, w, &wkopt, &lwork, &iwkopt,
                        &liwork, &info );
        lwork = (int)wkopt;
        work = (double*)malloc( lwork*sizeof(double) );
        liwork = iwkopt;
        iwork = (int*)malloc( liwork*sizeof(int) );

        clock_t start_dsyevd = clock();

        /* Solve eigenproblem */
        dsyevd_( "V", "L", &n, a, &lda, w, work, &lwork, iwork,
                        &liwork, &info );
        /* Check for convergence */
        if( info > 0 ) {
                printf( "The algorithm failed to compute eigenvalues.\n" );
                exit( 1 );
                free( (void*)a );
                free( (void*)iwork );
                free( (void*)work );
        }

        clock_t end_dsyevd = clock();
        const double time_dsyevd = static_cast<double>(end_dsyevd - start_dsyevd) / CLOCKS_PER_SEC * 1000.0;
        printf("dsyevd time %lf[ms]\n", time_dsyevd);

        ofstream time("result_run_time.txt", ios::app);    //txtファイル書き出し
        time << "dsyevd time " << time_dsyevd << "[ms]\n";
        time.close();

        clock_t start_kakiutusi = clock();
        Eigen::VectorXd eig_val = Eigen::VectorXd::Zero(N);
        Eigen::MatrixXd eig_vec = Eigen::MatrixXd::Zero(N, N);
        for (int ii = 0; ii < N; ii++) {
            eig_val(ii) = w[ii];
        }
    
        for(int ii = 0; ii < N; ii++ ){
            for(int jj = 0; jj < N; jj++ ){
                eig_vec(ii, jj) = a[ii+jj*N];
            }
        }
        
        clock_t end_kakiutusi = clock();
        const double time_kakiutusi = static_cast<double>(end_kakiutusi - start_kakiutusi) / CLOCKS_PER_SEC * 1000.0;
        printf("kakiutusi time %lf[ms]\n", time_kakiutusi);

        /* Print eigenvalues */
        //print_matrix( "Eigenvalues", 1, n, w, 1 );
        /* Print eigenvectors */
        //print_matrix( "Eigenvectors (stored columnwise)", n, n, a, lda );
        /* Free workspace */
        free( (void*)iwork );
        free( (void*)work );
        //exit( 0 );

    return forward_as_tuple(eig_val, eig_vec);
} /* End of DSYEVD Example */
