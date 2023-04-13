#pragma once
#include <stdlib.h>
#include <stdio.h>
#include "Parameter.h"
#define PRINT_MAT(X) cout << #X << ":\n" << X << endl << endl   //行列の出力

/* Complex datatype */
struct _dcomplex { double re, im; };
typedef struct _dcomplex dcomplex;

/* ZHEEVR prototype */
extern "C" {
   void zheevr_( const char& jobz, const char& range, const char& uplo, int* n, dcomplex* a,
                int* lda, double* vl, double* vu, int* il, int* iu, double* abstol,
                int* m, double* w, dcomplex* z, int* ldz, int* isuppz, dcomplex* work,
                int* lwork, double* rwork, int* lrwork, int* iwork, int* liwork,
                int* info );
};
/* Auxiliary routines prototypes */
//extern void print_matrix( const char* desc, int m, int n, dcomplex* a, int lda );
//extern void print_rmatrix( const char* desc, int m, int n, double* a, int lda );

/* Auxiliary routine: printing a matrix */
void print_matrix( const char* desc, int m, int n, dcomplex* a, int lda ) {
        int i, j;
        printf( "\n %s\n", desc );
        for( i = 0; i < m; i++ ) {
                for( j = 0; j < n; j++ )
                        printf( "\t%6.3f", a[i+j*lda].re );
                printf( "\n" );
        }
}

/* Auxiliary routine: printing a real matrix */
void print_rmatrix( const char* desc, int m, int n, double* a, int lda ) {
        int i, j;
        printf( "\n %s\n", desc );
        for( i = 0; i < m; i++ ) {
                for( j = 0; j < n; j++ ) printf( " %6.6f", a[i+j*lda] );
                printf( "\n" );
        }
}

/* Main program */
tuple<Eigen::VectorXd, Eigen::Matrix<complex<double>, Eigen::Dynamic, Eigen::Dynamic>> diagonalization_zheevr(Eigen::Matrix<complex<double>, Eigen::Dynamic, Eigen::Dynamic> mat) {
        /* Parameters */
        #define N mat.rows()
        #define LDA N
        #define LDZ N
        /* Locals */
        int n = N, lda = LDA, ldz = LDZ, il, iu, m, info, lwork, lrwork, liwork;
        double abstol, vl, vu;
        int iwkopt;
        int* iwork;
        double rwkopt;
        double* rwork;
        dcomplex wkopt;
        dcomplex* work;
        /* Local arrays */
        int* isuppz;
        double* w;
        dcomplex* z;
        dcomplex* a;
        isuppz = (int*)malloc( N*sizeof(int) );
        w = (double*)malloc( N*sizeof(double) );
        z = (dcomplex*)malloc( N*N*sizeof(dcomplex) );
        a = (dcomplex*)malloc( N*N*sizeof(dcomplex) );
        for(int i=0; i<N; i++ ){
            for(int j=0; j<=i; j++ ){
                a[i+j*lda].re  = mat(i, j).real();
                a[i+j*lda].im = mat(i, j).imag();
            }
        }

       
        /* Executable statements */
        /* Negative abstol means using the default value */
        abstol = -1.0;
        /* Set VL, VU to compute eigenvalues in half-open (VL,VU] interval */
        vl = -100.0;
        vu = 10000.0;
        /* Query and allocate the optimal workspace */
        lwork = -1;
        lrwork = -1;
        liwork = -1;
        zheevr_( 'V', 'V', 'L', &n, a, &lda, &vl, &vu, &il, &iu,
                        &abstol, &m, w, z, &ldz, isuppz, &wkopt, &lwork, &rwkopt, &lrwork,
                        &iwkopt, &liwork, &info );
        lwork = (int)wkopt.re;
        work = (dcomplex*)malloc( lwork*sizeof(dcomplex) );
        lrwork = (int)rwkopt;
        rwork = (double*)malloc( lrwork*sizeof(double) );
        liwork = iwkopt;
        iwork = (int*)malloc( liwork*sizeof(int) );

        clock_t start_zheevr = clock();

        /* Solve eigenproblem */
        zheevr_( 'V', 'V', 'L', &n, a, &lda, &vl, &vu, &il, &iu,
                        &abstol, &m, w, z, &ldz, isuppz, work, &lwork, rwork, &lrwork,
                        iwork, &liwork, &info );
        /* Check for convergence */
        if( info > 0 ) {
                printf( "The algorithm failed to compute eigenvalues.\n" );
                free( (void*)isuppz );
                free( (void*)w );
                free( (void*)z );
                free( (void*)a );
                free( (void*)iwork );
                free( (void*)rwork );
                free( (void*)work );
                exit( 1 );

        }

        clock_t end_zheevr = clock();
        const double time_zheevr = static_cast<double>(end_zheevr - start_zheevr) / CLOCKS_PER_SEC;
        printf("zheevr time %lf[s]\n", time_zheevr);

        ofstream time("result_run_time.txt", ios::app);    //txtファイル書き出し
        time << "zheevr time " << time_zheevr << "[s]\n";
        time.close();

        Eigen::VectorXd eig_val = Eigen::VectorXd::Zero(N);
        Eigen::Matrix<complex<double>, Eigen::Dynamic, Eigen::Dynamic> eig_vec = Eigen::Matrix<complex<double>, Eigen::Dynamic, Eigen::Dynamic>::Zero(N, N);
        for (int ii = 0; ii < N; ii++) {
            eig_val(ii) = w[ii];
        }
    
        for(int ii = 0; ii < N; ii++ ){
            for(int jj = 0; jj < N; jj++ ){
                eig_vec(ii, jj).real(z[ii+jj*N].re);
                eig_vec(ii, jj).imag(z[ii+jj*N].im);
            }
        }
        
        /* Print the number of eigenvalues found */
        //printf( "\n The total number of eigenvalues found:%2i\n", m );
        /* Print eigenvalues */
        //print_rmatrix( "Selected eigenvalues", 1, m, w, 1 );
        /* Print eigenvectors */
        //print_matrix( "Selected eigenvectors (stored columnwise)", n, m, z, ldz );
        /* Free workspace */
        free( (void*)isuppz );
        free( (void*)w );
        free( (void*)z );
        free( (void*)a );
        free( (void*)iwork );
        free( (void*)rwork );
        free( (void*)work );
        //exit( 0 );
        

    return forward_as_tuple(eig_val, eig_vec);
} /* End of ZHEEVR Example */
