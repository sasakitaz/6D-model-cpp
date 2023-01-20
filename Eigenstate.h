#pragma once
#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>  // abs() for integer
#include <numeric>  //accumulate
#include <cmath>
#include <bits/stdc++.h>    //pi: M_PI 
#include "Parameter.h"  
#define PRINT_MAT(X) cout << #X << ":\n" << X << endl << endl 

using namespace std;


//for non symmetry matrix (real)
int eigen_state(int m, Eigen::MatrixXd eigen_value, Eigen::MatrixXd eigen_vector) 
{

    char fname[30];
    sprintf(fname, "result_value_%d.txt", m);
    ofstream e_value(fname);
    e_value << "n0 = " << n0 << "\t" << "npls = " << npls << "\t" << "nmns = " << nmns << "\t" << "j = " << j << "\t" << "m = " << calcm  << "\n" ;
    e_value << "internal rotation parameter, (V3, V4) = " << "(" << V3 << ", "  << V4 << ") \n";
    e_value << "stretch and bend parameter, (kzz, az, kxx, kxxz) = " << "(" << kzz << ", "  << az << ", "  << kxx << ", "  << kxxz << ") \n";
    e_value << "R coupling parameter, (VR3, VRR3, VR4, VRR4) = " << "(" << VR3 << ", "  << VRR3 << ", "  << VR4 << ", "  << VRR4 << ") \n";
    e_value << "rho coupling parameter, (Vrho3, Vrho4) = " << "(" << Vrho3 << ", "  << Vrho4 << ") \n";
    e_value << "angular momentum coupling, (Vam31, VamR31) = " << "(" << Vam31 << ", "  << VamR31 << ") \n\n";
    e_value << eigen_value.transpose();
    e_value.close();
    
    char ffname[30];
    sprintf(ffname, "result_vector_%d.txt", m);
    ofstream e_vector(ffname); 
    e_vector << "n0 = " << n0 << "\t" << "npls = " << npls << "\t" << "nmns = " << nmns << "\t" << "j = " << j << "\t" << "m = " << calcm  << "\n" ;
    e_vector << "internal rotation parameter, (V3, V4) = " << "(" << V3 << ", " << V4 << ") \n";
    e_vector << "stretch and bend parameter, (kzz, az, kxx, kxxz) = " << "(" << kzz << ", "  << az << ", "  << kxx << ", "  << kxxz << ") \n";
    e_vector << "R coupling parameter, (VR3, VRR3, VR4, VRR4) = " << "(" << VR3 << ", "  << VRR3 << ", "  << VR4 << ", "  << VRR4 << ") \n";
    e_vector << "rho coupling parameter, (Vrho3, Vrho4) = " << "(" << Vrho3 << ", "  << Vrho4 << ") \n";
    e_vector << "angular momentum coupling, (Vam31, VamR31) = " << "(" << Vam31 << ", "  << VamR31 << ") \n\n";

    eigen_vector.transposeInPlace();    //transpose unitary matrix
    int jjcount = 0;
    vector<int> j2count;
    for (int ja = 0;ja <= j; ja++) {
        jjcount += 2*ja + 1;
        j2count.push_back(jjcount);
    }
    int quantn0;
    int quantv;
    int quantl;
    int quantj;
    int quantk;
    int quantmB;
    int quantm;
    int csurp;
    int jsurp;
    int ksurp;
    vector<vector<int>> quant;  //debag
    vector<int> dim;
    vector<int> allowl;
    vector<int> allowmB;
    MatrixParameter hoge(m, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
    int dimension = hoge.matrix_size();    //dimension number
    
    for (int mB = - j; mB <= j; mB++) {
        for (int ll = - nmns; ll <= npls; ll++) {
            int jdim = 0;
            for (int num = 0; num <= j; num++) {
                jdim += 2*abs(num) + 1;
            }
            int mminus = 0;
            for (int mm = 0; mm <= abs(mB); mm++){
                if (mm == 0){mminus += 0;}
                else {mminus += 2*(abs(mm) - 1) + 1;}
            }
            if (ll + mB == m) {
                dim.push_back((npls + 1 - abs(ll))*(jdim - mminus));
                allowl.push_back(ll);
                allowmB.push_back(mB);
                
            }
        }
    }
    int sum_dim = accumulate(dim.begin(), dim.end(), 0);
    auto size_dim = dim.size(); 

    for (int r = 0; r <= eigen_vector.rows() - 1; r++) {
        double cdouble = 0;
        for (int c = 0; c <= eigen_vector.cols() - 1; c++) {
            if (abs(eigen_vector(r, c)) > diag_coeff) {
                quantn0 = c / sum_dim;  //quantum number n0
                quantm = m; //quantum number m
                csurp = c % sum_dim; 

                vector<int> sumd(1,0);
                //quantum number npls, nmns, mB
                for (int d = 0; d <= size_dim - 1; d++) {
                    sumd.push_back(sumd[d] + dim[d]);
                    if (sumd[d] <= csurp && csurp < sumd[d + 1]) {
                        quantl = allowl[d];
                        quantmB = allowmB[d];
                        jsurp = csurp - sumd[d];    //jsurp: number of each (mB, l) matrix

                        if (quantmB == 0) {
                            quantv = 2*(jsurp / jjcount) + abs(quantl);
                            ksurp = jsurp % j2count[j];
                        }
                        else if ( j == 0) {
                            quantv = 2*(jsurp / jjcount) + abs(quantl);
                            ksurp = jsurp % j2count[j];
                        }
                        else {
                            quantv = 2*(jsurp / (jjcount - j2count[abs(quantmB) - 1])) + abs(quantl);
                            ksurp = jsurp % (j2count[j] - j2count[abs(quantmB) - 1]) + j2count[abs(quantmB) - 1];
                        }

                        if (ksurp == 0) {
                            quantj = abs(quantmB);
                            quantk = quantmB;

                            cdouble += abs(eigen_vector(r, c))*abs(eigen_vector(r, c)); //debag
                            
                            if (eigen_value(r) < 300 && quantn0 < 5 && quantv < 4 && abs(quantl) < 4 && quantj < 6 && abs(quantk) < 6 && abs(quantl) < 6){

                            e_vector << eigen_value(r) << ",\t" 
                                     << eigen_vector(r, c) << ",\t" 
                                     << quantn0 << ",\t" 
                                     << quantv << ",\t" 
                                     << quantl << ",\t" 
                                     << quantj << ",\t" 
                                     << quantk << ",\t" 
                                     << quantmB << ",\t" 
                                     << quantm << "\n" ;
                            }
                        }
                        for (int jc = 0; jc <= j; jc++) {
                            if (ksurp >= j2count[jc] && ksurp < j2count[jc + 1]) {
                                quantj = jc + 1;
                                quantk = ksurp - j2count[jc] - quantj;
                                
                                cdouble += abs(eigen_vector(r, c))*abs(eigen_vector(r, c)); //debag

                                if (eigen_value(r) < 300 && quantn0 < 5 && quantv < 4 && abs(quantl) < 4 && quantj < 6 && abs(quantk) < 6 && abs(quantl) < 6){

                                e_vector << eigen_value(r) << ",\t" 
                                         << eigen_vector(r, c) << ",\t" 
                                         << quantn0 << ",\t" 
                                         << quantv << ",\t" 
                                         << quantl << ",\t" 
                                         << quantj << ",\t" 
                                         << quantk << ",\t" 
                                         << quantmB << ",\t" 
                                         << quantm << "\n" ;
                                }
                            }
                        }
                    }
                }

            }
        }
    }
    e_vector.close();
    return 0;
}

//For symmetrized matrix (complex)
int eigen_state(int m, Eigen::MatrixXd eigen_value, Eigen::Matrix<complex<double>, Eigen::Dynamic, Eigen::Dynamic> eigen_vector, Eigen::Matrix<complex<double>, Eigen::Dynamic, Eigen::Dynamic> U, string sym) 
{
    ofstream e_value("result_" + sym + "value.txt"); 
    e_value << "n0 = " << n0 << "\t" << "npls = " << npls << "\t" << "nmns = " << nmns << "\t" << "j = " << j << "\t" << "m = " << calcm  << "\n" ;
    e_value << "internal rotation parameter, (V3, V4) = " << "(" << V3 << ", "  << V4 << ") \n";
    e_value << "stretch and bend parameter, (kzz, az, kxx, kxxz) = " << "(" << kzz << ", "  << az << ", "  << kxx << ", "  << kxxz << ") \n";
    e_value << "R coupling parameter, (VR3, VRR3, VR4, VRR4) = " << "(" << VR3 << ", "  << VRR3 << ", "  << VR4 << ", "  << VRR4 << ") \n";
    e_value << "rho coupling parameter, (Vrho3, Vrho4) = " << "(" << Vrho3 << ", "  << Vrho4 << ") \n";
    e_value << "angular momentum coupling, (Vam31, VamR31) = " << "(" << Vam31 << ", "  << VamR31 << ") \n\n";
    e_value << eigen_value.transpose();
    e_value.close();
    
    ofstream e_vector("result_" + sym + "vector.txt");
    e_vector << "n0 = " << n0 << "\t" << "npls = " << npls << "\t" << "nmns = " << nmns << "\t" << "j = " << j << "\t" << "m = " << calcm  << "\n" ;
    e_vector << "internal rotation parameter, (V3, V4) = " << "(" << V3 << ", " << V4 << ") \n";
    e_vector << "stretch and bend parameter, (kzz, az, kxx, kxxz) = " << "(" << kzz << ", "  << az << ", "  << kxx << ", "  << kxxz << ") \n";
    e_vector << "R coupling parameter, (VR3, VRR3, VR4, VRR4) = " << "(" << VR3 << ", "  << VRR3 << ", "  << VR4 << ", "  << VRR4 << ") \n";
    e_vector << "rho coupling parameter, (Vrho3, Vrho4) = " << "(" << Vrho3 << ", "  << Vrho4 << ") \n";
    e_vector << "angular momentum coupling, (Vam31, VamR31) = " << "(" << Vam31 << ", "  << VamR31 << ") \n\n";

    eigen_vector = U.conjugate().transpose()*eigen_vector;
    eigen_vector.transposeInPlace();   //transpose unitary matrix
    int jjcount = 0;
    vector<int> j2count;
    for (int ja = 0;ja <= j; ja++) {
        jjcount += 2*ja + 1;
        j2count.push_back(jjcount);
    }
    int quantn0;
    int quantv;
    int quantl;
    int quantj;
    int quantk;
    int quantmB;
    int quantm;
    int csurp;
    int jsurp;
    int ksurp;
    vector<vector<int>> quant;  //debag
    vector<int> dim;
    vector<int> allowl;
    vector<int> allowmB;
    MatrixParameter hoge(m, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
    int dimension = hoge.matrix_size();    //dimension number
    
    for (int mB = - j; mB <= j; mB++) {
        for (int ll = - nmns; ll <= npls; ll++) {
            int jdim = 0;
            for (int num = 0; num <= j; num++) {
                jdim += 2*abs(num) + 1;
            }
            int mminus = 0;
            for (int mm = 0; mm <= abs(mB); mm++){
                if (mm == 0){mminus += 0;}
                else {mminus += 2*(abs(mm) - 1) + 1;}
            }
            if (ll + mB == m) {
                dim.push_back((npls + 1 - abs(ll))*(jdim - mminus));
                allowl.push_back(ll);
                allowmB.push_back(mB);
                
            }
        }
    }
    int sum_dim = accumulate(dim.begin(), dim.end(), 0);
    auto size_dim = dim.size();

    for (int r = 0; r <= eigen_vector.rows() - 1; r++) {
        double cdouble = 0;
        for (int c = 0; c <= eigen_vector.cols() - 1; c++) {
            if (abs(eigen_vector(r, c)) > 0.2) {
                quantn0 = c / sum_dim;  //quantum number n0
                quantm = m; //quantum number m
                csurp = c % sum_dim;

                vector<int> sumd(1,0);
                //quantum number npls, nmns, mB
                for (int d = 0; d <= size_dim - 1; d++) {
                    sumd.push_back(sumd[d] + dim[d]);

                    if (sumd[d] <= csurp && csurp < sumd[d + 1]) {
                        quantl = allowl[d];
                        quantmB = allowmB[d];
                        jsurp = csurp - sumd[d];    //jsurp: number of each (mB, l) matrix

                        if (quantmB == 0) {
                            quantv = 2*(jsurp / jjcount) + abs(quantl);
                            ksurp = jsurp % j2count[j];
                        }
                        else if ( j == 0) {
                            quantv = 2*(jsurp / jjcount) + abs(quantl);
                            ksurp = jsurp % j2count[j];
                        }
                        else {
                            quantv = 2*(jsurp / (jjcount - j2count[abs(quantmB) - 1])) + abs(quantl);
                            ksurp = jsurp % (j2count[j] - j2count[abs(quantmB) - 1]) + j2count[abs(quantmB) - 1];
                        }

                        if (ksurp == 0) {
                            quantj = abs(quantmB);
                            quantk = quantmB;

                            cdouble += abs(eigen_vector(r, c))*abs(eigen_vector(r, c)); //debag
                            
                            e_vector << eigen_value(r) << ",\t" 
                                     << eigen_vector(r, c) << ",\t"
                                     << abs(eigen_vector(r, c)) << ",\t" 
                                     << quantn0 << ",\t" 
                                     << quantv << ",\t" 
                                     << quantl << ",\t" 
                                     << quantj << ",\t" 
                                     << quantk << ",\t" 
                                     << quantmB << ",\t" 
                                     << quantm << "\n" ;
                        }
                        for (int jc = 0; jc <= j; jc++) {
                            if (ksurp >= j2count[jc] && ksurp < j2count[jc + 1]) {

                                quantj = jc + 1;
                                quantk = ksurp - j2count[jc] - quantj;
                                
                                cdouble += abs(eigen_vector(r, c))*abs(eigen_vector(r, c)); //debag

                                e_vector << eigen_value(r) << ",\t" 
                                         << eigen_vector(r, c) << ",\t"
                                         << abs(eigen_vector(r, c)) << ",\t" 
                                         << quantn0 << ",\t" 
                                         << quantv << ",\t" 
                                         << quantl << ",\t" 
                                         << quantj << ",\t" 
                                         << quantk << ",\t" 
                                         << quantmB << ",\t" 
                                         << quantm << "\n" ;
                            }
                        }
                    }
                }

            }
        }
    }
    e_vector.close();
    return 0;
}