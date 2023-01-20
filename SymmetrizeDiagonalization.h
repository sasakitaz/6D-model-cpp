#pragma once
#include <iostream>
#include <vector>
#include <cstdlib>  // abs() for integer
#include <numeric>  //accumulate
#include <cmath>
#include <time.h>
#include "Symmetrize.h"
#include "MatrixElement.h"
#include "Eigenstate.h" 
#include "diagonalization_zheevr.h" 
#include "diagonalization_dsyevd.h"

#define PRINT_MAT(X) cout << #X << ":\n" << X << endl << endl   //行列の出力

int symmetrized_diagonalization6(int m, Eigen::MatrixXd H) {
    int q = 0;
    vector<string> project_name{"A1", "E1", "E2", "F1"};
    MatrixParameter dim(m, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
    vector<complex<double>> zero(dim.matrix_size(), 0);
    Eigen::Vector<vector<vector<complex<double>>>, 4> projects; //from Symmetrize.h
    projects = Projects6(m);
    for (auto& project : projects) {
        q += 1;
        vector<vector<complex<double>>> project_nonzero = project;

        int dim_nonzero = project.size();
        cout << "dimension, " << dim_nonzero << endl;
        if (dim_nonzero != 0) {

            //vector -> Eigen
            clock_t start_rewrite = clock();
            Eigen::Matrix<complex<double>, Eigen::Dynamic, Eigen::Dynamic> project_Eigen = Eigen::Matrix<complex<double>, Eigen::Dynamic, Eigen::Dynamic>::Zero(project_nonzero.size(), project_nonzero[0].size());
            for (int row = 0; row < project_nonzero.size(); row++){
                for(int column = 0; column < project_nonzero[0].size(); column++){
                    project_Eigen(row, column) = project[row][column];
                }
            }
            clock_t end_rewrite = clock();
            double time_rewrite = double(end_rewrite - start_rewrite) / CLOCKS_PER_SEC;
            ofstream time("result_run_time.txt", ios::app);    //txtファイル書き出し
            if      (q == 1){printf("A1 rewrite time %lf[s]\n", time_rewrite); time << "A1 rewrite time " << time_rewrite << "[s]\n";}
            else if (q == 2){printf("E1 rewrite time %lf[s]\n", time_rewrite); time << "E1 rewrite time " << time_rewrite << "[s]\n";}
            else if (q == 3){printf("E2 rewrite time %lf[s]\n", time_rewrite); time << "E2 rewrite time " << time_rewrite << "[s]\n";}
            else if (q == 4){printf("F1 rewrite time %lf[s]\n", time_rewrite); time << "F1 rewrite time " << time_rewrite << "[s]\n";}
            time.close();

            Eigen::Matrix<complex<double>, Eigen::Dynamic, Eigen::Dynamic> H_Gamma;
            H_Gamma = project_Eigen*H*project_Eigen.conjugate().transpose();
            if (dim_nonzero != 0){
                clock_t start_symdiag = clock();
                //diagonalization
                Eigen::VectorXd eig_val = Eigen::VectorXd::Zero(H_Gamma.rows());
                Eigen::Matrix<complex<double>, Eigen::Dynamic, Eigen::Dynamic> eig_vec = Eigen::Matrix<complex<double>, Eigen::Dynamic, Eigen::Dynamic>::Zero(H_Gamma.rows(), H_Gamma.cols());
                tie(eig_val, eig_vec) = diagonalization_zheevr(H_Gamma);    //diagonalization, complex Hermite matrix from diagonalization_zheevr.h
                eigen_state(m, eig_val, eig_vec, project_Eigen, project_name[q - 1]); //calculation of eigenvector's component from Eigenstate.h
                clock_t end_symdiag = clock();
                double time_symdiag = double(end_symdiag - start_symdiag) / CLOCKS_PER_SEC;
                ofstream time("result_run_time.txt", ios::app);
                if      (q == 1){printf("A1 symmetrized diagonalizatoin time %lf[s]\n", time_symdiag); time << "A1 symmetrized diagonalizatoin time " << time_symdiag << "[s]\n";}
                else if (q == 2){printf("E1 symmetrized diagonalizatoin time %lf[s]\n", time_symdiag); time << "E1 symmetrized diagonalizatoin time " << time_symdiag << "[s]\n";}
                else if (q == 3){printf("F1 symmetrized diagonalizatoin time %lf[s]\n", time_symdiag); time << "F1 symmetrized diagonalizatoin time " << time_symdiag << "[s]\n";}
                time.close();
            } 
        }
        else{
            vector<vector<complex<double>>> project(dim.matrix_size(), vector<complex<double>>(dim.matrix_size(), 0)); //diagonarization from Eigen
            for (int i = 0; i < dim.matrix_size(); i++) {project[i][i] = 1;}
        }
    }
    return 0;
}

int symmetrized_diagonalization3(int m, Eigen::MatrixXd H) {
    int q = 0;
    vector<string> project_name{"A", "E", "F"};
    MatrixParameter dim(m, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
    vector<complex<double>> zero(dim.matrix_size(), 0);
    Eigen::Vector<vector<vector<complex<double>>>, 3> projects; //from Symmetrize.h
    tie(ignore, projects, ignore) = Projects(m);
    for (auto& project : projects) {
        q += 1;
        cout << "dimension, " << project.size() << endl;
        if (project.size() != 0) {
            //vector -> Eigen
            clock_t start_rewrite = clock();
            Eigen::Matrix<complex<double>, Eigen::Dynamic, Eigen::Dynamic> project_Eigen = Eigen::Matrix<complex<double>, Eigen::Dynamic, Eigen::Dynamic>::Zero(project.size(), project[0].size());
            for (int row = 0; row < project.size(); row++){
                for(int column = 0; column < project[0].size(); column++){
                    project_Eigen(row, column) = project[row][column];
                }
            }
            clock_t end_rewrite = clock();
            double time_rewrite = double(end_rewrite - start_rewrite) / CLOCKS_PER_SEC;

            ofstream time("result_run_time.txt", ios::app);
            if      (q == 1){printf("A rewrite time %lf[s]\n", time_rewrite); time << "A rewrite time " << time_rewrite << "[s]\n";}
            else if (q == 2){printf("E rewrite time %lf[s]\n", time_rewrite); time << "E rewrite time " << time_rewrite << "[s]\n";}
            else if (q == 3){printf("F rewrite time %lf[s]\n", time_rewrite); time << "F rewrite time " << time_rewrite << "[s]\n";}
            time.close();

            Eigen::Matrix<complex<double>, Eigen::Dynamic, Eigen::Dynamic> H_Gamma;
            H_Gamma = project_Eigen*H*project_Eigen.conjugate().transpose();
            if (project.size() != 0){
                clock_t start_symdiag = clock();
                //diagonalization
                Eigen::VectorXd eig_val = Eigen::VectorXd::Zero(H_Gamma.rows());
                Eigen::Matrix<complex<double>, Eigen::Dynamic, Eigen::Dynamic> eig_vec = Eigen::Matrix<complex<double>, Eigen::Dynamic, Eigen::Dynamic>::Zero(H_Gamma.rows(), H_Gamma.cols());
                tie(eig_val, eig_vec) = diagonalization_zheevr(H_Gamma);    //diagonalization, complex Hermite matrix from diagonalization_zheevr.h
                eigen_state(m, eig_val, eig_vec, project_Eigen, project_name[q - 1]); //calculation of eigenvector's component from Eigenstate.h
                clock_t end_symdiag = clock();
                double time_symdiag = double(end_symdiag - start_symdiag) / CLOCKS_PER_SEC;
                ofstream time("result_run_time.txt", ios::app);    //txtファイル書き出し
                if      (q == 1){printf("A symmetrized diagonalizatoin time %lf[s]\n", time_symdiag); time << "A symmetrized diagonalizatoin time " << time_symdiag << "[s]\n";}
                else if (q == 2){printf("E symmetrized diagonalizatoin time %lf[s]\n", time_symdiag); time << "E symmetrized diagonalizatoin time " << time_symdiag << "[s]\n";}
                else if (q == 3){printf("F symmetrized diagonalizatoin time %lf[s]\n", time_symdiag); time << "F symmetrized diagonalizatoin time " << time_symdiag << "[s]\n";}
                time.close();
            } 
        }
        else{
            vector<vector<complex<double>>> project(dim.matrix_size(), vector<complex<double>>(dim.matrix_size(), 0)); //diagonalization from Eigen
            for (int i = 0; i < dim.matrix_size(); i++) {project[i][i] = 1;}
        }
    }
    return 0;
}

int symmetrized_diagonalization1(int m, Eigen::MatrixXd H) {
    int q = 0;
    MatrixParameter dim(m, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
    vector<complex<double>> zero(dim.matrix_size(), 0);
    Eigen::Vector<vector<vector<complex<double>>>, 1> projects; //from Symmetrize.h
    tie(ignore, ignore, projects) = Projects(m);
    vector<vector<complex<double>>> project = projects(0);
    int dim_nonzero = project.size();
    //vector -> Eigen
    Eigen::Matrix<complex<double>, Eigen::Dynamic, Eigen::Dynamic> project_Eigen = Eigen::Matrix<complex<double>, Eigen::Dynamic, Eigen::Dynamic>::Zero(project.size(), project[0].size());
    for (int row = 0; row < project.size(); row++){
        for(int column = 0; column < project[0].size(); column++){
            project_Eigen(row, column) = project[row][column];
        }
    }
    
    Eigen::Matrix<complex<double>, Eigen::Dynamic, Eigen::Dynamic> H_Gamma;
    H_Gamma = project_Eigen*H*project_Eigen.conjugate().transpose();
    H_Gamma = H_Gamma(Eigen::seqN(0, dim_nonzero), Eigen::seqN(0, dim_nonzero));
    if (dim_nonzero != 0){
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix<complex<double>, Eigen::Dynamic, Eigen::Dynamic>> ev(H_Gamma); //diagonalization from Eigen
        eigen_state(m, ev.eigenvalues(), ev.eigenvectors(), project_Eigen, ""); //calculation of eigenvector's component from Eigenstate.h
    }
    return 0;
}

int diagonalization(int m, Eigen::MatrixXd H) {
    //diagonalization
    Eigen::VectorXd eig_val = Eigen::VectorXd::Zero(H.rows());
    Eigen::MatrixXd eig_vec = Eigen::MatrixXd::Zero(H.rows(), H.rows());
    tie(eig_val, eig_vec) = diagonalization_dsyevd(H);    //diagonalization, real symmetry matrix from diagonalization_dsyevd.h
    eigen_state(m, eig_val, eig_vec);     //calculation of eigenvector's component from Eigenstate.h
    return 0;
}
