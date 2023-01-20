#include <iostream>
#include <vector>
#include <cstdlib>  // abs() for integer
#include <numeric>  //accumulate
#include <cmath>
#include <complex.h>
#include <algorithm>
#include <bits/stdc++.h>    //pi: M_PI
#include "MatrixElement.h"
#include "WignerdSmall_pihalf.h"   //Wigner d function

#define PRINT_MAT(X) cout << #X << ":\n" << X << endl << endl   //行列の出力

using namespace std;

class SymmetrizeMatrix : public MatrixElement
{
public:
public:
    SymmetrizeMatrix(int temm, int temjrow, int temmBrow, int temkrow, int temn0row, int temvrow, int temjcolumn, int temmBcolumn, int temkcolumn, int temn0column, int temvcolumn) : MatrixElement(temm, temjrow, temmBrow, temkrow, temn0row, temvrow, temjcolumn, temmBcolumn, temkcolumn, temn0column, temvcolumn)
    {
    }
    double A1();
    double E1();
    double E2();
    double F1();
    complex<double> F2();
    complex<double> F3();

};

double SymmetrizeMatrix::A1(){
    double me_A1;
    if (krow%2 == 0) {
        if (jrow == jcolumn && n0row == n0column && nplsrow == nplscolumn && nmnsrow == nmnsrow && mBrow == mBcolumn) {
            if (krow == 0 && kcolumn == 0) {
                me_A1 = + 1
                        + pow((-1), jrow)*pow((-1), (krow/2))
                        + 2*wigner_d_small_pihalf(jrow)(jrow + kcolumn, jrow + krow)*cos(M_PI*(kcolumn + krow)/4)
                        + 2*wigner_d_small_pihalf(jrow)(jrow - kcolumn, jrow + krow)*cos(M_PI*(kcolumn + krow)/4)*pow((-1), jrow)*pow((-1), (kcolumn/2));
                if (abs(me_A1) < pow(10, -10)) { me_A1 = 0; }
                        
            }
            else if (krow == kcolumn) {
                me_A1 = + 1
                        + 2*wigner_d_small_pihalf(jrow)(jrow + kcolumn, jrow + krow)*cos(M_PI*(kcolumn + krow)/4)
                        + 2*wigner_d_small_pihalf(jrow)(jrow - kcolumn, jrow + krow)*cos(M_PI*(- kcolumn + krow)/4)*pow((-1), jrow)*pow((-1), (- kcolumn/2));
                if (abs(me_A1) < pow(10, -10)) { me_A1 = 0; }
            }
            else if (krow == - kcolumn) {
                me_A1 = + pow((-1), jrow)*pow((-1), (krow/2))
                        + 2*wigner_d_small_pihalf(jrow)(jrow + kcolumn, jrow + krow)*cos(M_PI*(kcolumn + krow)/4)
                        + 2*wigner_d_small_pihalf(jrow)(jrow - kcolumn, jrow + krow)*cos(M_PI*(- kcolumn + krow)/4)*pow((-1), jrow)*pow((-1), (- kcolumn/2));
                if (abs(me_A1) < pow(10, -10)) { me_A1 = 0; }
            }
            else if (kcolumn%2 == 0) {
                me_A1 = + 2*wigner_d_small_pihalf(jrow)(jrow + kcolumn, jrow + krow)*cos(M_PI*(kcolumn + krow)/4)
                        + 2*wigner_d_small_pihalf(jrow)(jrow - kcolumn, jrow + krow)*cos(M_PI*(- kcolumn + krow)/4)*pow((-1), jrow)*pow((-1), (- kcolumn/2));
                if (abs(me_A1) < pow(10, -10)) { me_A1 = 0; }
            }
            else {
                me_A1 = 0;
            }
        }
        else {
            me_A1 = 0;
        }
    }
    else {
        me_A1 = 0;
    }
    return me_A1;
}

double SymmetrizeMatrix::E1() {
    double me_E1;
    if (krow%2 == 0) {
        if (jrow == jcolumn && n0row == n0column && nplsrow == nplscolumn && nmnsrow == nmnsrow && mBrow == mBcolumn) {
            if (krow == 0 && kcolumn == 0) {
                me_E1 = + 1 
                        + pow((-1), jrow)*pow((-1), (krow/2)) 
                        + 2*wigner_d_small_pihalf(jrow)(jrow + kcolumn, jrow + krow)*cos(((kcolumn + krow)*M_PI/4 - 2*M_PI/3))
                        + 2*wigner_d_small_pihalf(jrow)(jrow - kcolumn, jrow + krow)*cos(((kcolumn + krow)/4 - 2*M_PI/3))*pow((-1), jrow)*pow((-1), (kcolumn/2));
                if (abs(me_E1) < pow(10, -10)) { me_E1 = 0; }
            }
            else if (krow == kcolumn) {
                me_E1 = + 1
                        + 2*wigner_d_small_pihalf(jrow)(jrow + kcolumn, jrow + krow)*cos(((kcolumn + krow)*M_PI/4 - 2*M_PI/3))
                        + 2*wigner_d_small_pihalf(jrow)(jrow - kcolumn, jrow + krow)*cos(((- kcolumn + krow)*M_PI/4 - 2*M_PI/3))*pow((-1), jrow)*pow((-1), (- kcolumn/2));
                if (abs(me_E1) < pow(10, -10)) { me_E1 = 0; }
            }
            else if (krow == - kcolumn) {
                me_E1 = + pow((-1), jrow)*pow((-1), (krow/2))
                        + 2*wigner_d_small_pihalf(jrow)(jrow + kcolumn, jrow + krow)*cos(((kcolumn + krow)*M_PI/4 - 2*M_PI/3))
                        + 2*wigner_d_small_pihalf(jrow)(jrow - kcolumn, jrow + krow)*cos(((- kcolumn + krow)*M_PI/4 - 2*M_PI/3))*pow((-1), jrow)*pow((-1), (- kcolumn/2));
                if (abs(me_E1) < pow(10, -10)) { me_E1 = 0; }                                                                                        
            }
            else if (kcolumn%2 == 0) {
                me_E1 = + 2*wigner_d_small_pihalf(jrow)(jrow + kcolumn, jrow + krow)*cos(((kcolumn + krow)*M_PI/4 - 2*M_PI/3))
                        + 2*wigner_d_small_pihalf(jrow)(jrow - kcolumn, jrow + krow)*cos(((- kcolumn + krow)*M_PI/4 - 2*M_PI/3))*pow((-1), jrow)*pow((-1), (- kcolumn/2));
                if (abs(me_E1) < pow(10, -10)) { me_E1 = 0; }
            }
            else {
                me_E1 = 0;
            }
        }
        else {
            me_E1 = 0;
        }
    }
    else {
        me_E1 = 0;
    }
    return me_E1;
}

double SymmetrizeMatrix::E2() {
    double me_E2;
    if (krow%2 == 0) {
        if (jrow == jcolumn && n0row == n0column && nplsrow == nplscolumn && nmnsrow == nmnsrow && mBrow == mBcolumn) {
            if (krow == 0 && kcolumn == 0) {
                me_E2 = + 1
                        + pow((-1), jrow)*pow((-1), (krow/2))
                        + 2*wigner_d_small_pihalf(jrow)(jrow + kcolumn, jrow + krow)*cos(((kcolumn + krow)*M_PI/4 + 2*M_PI/3))
                        + 2*wigner_d_small_pihalf(jrow)(jrow - kcolumn, jrow + krow)*cos(((kcolumn + krow)*M_PI/4 + 2*M_PI/3))*pow((-1), jrow)*pow((-1), (kcolumn/2));
                if (abs(me_E2) < pow(10, -10)) { me_E2 = 0; }
            }
            else if (krow == kcolumn) {
                me_E2 = + 1
                        + 2*wigner_d_small_pihalf(jrow)(jrow + kcolumn, jrow + krow)*cos(((kcolumn + krow)*M_PI/4 + 2*M_PI/3))
                        + 2*wigner_d_small_pihalf(jrow)(jrow - kcolumn, jrow + krow)*cos(((- kcolumn + krow)*M_PI/4 + 2*M_PI/3))*pow((-1), jrow)*pow((-1), (- kcolumn/2));
                if (abs(me_E2) < pow(10, -10)) { me_E2 = 0; }
            }
            else if (krow == - kcolumn) {
                me_E2 = + pow((-1), jrow)*pow((-1), (krow/2))
                        + 2*wigner_d_small_pihalf(jrow)(jrow + kcolumn, jrow + krow)*cos(((kcolumn + krow)*M_PI/4 + 2*M_PI/3))
                        + 2*wigner_d_small_pihalf(jrow)(jrow - kcolumn, jrow + krow)*cos(((- kcolumn + krow)*M_PI/4 + 2*M_PI/3))*pow((-1), jrow)*pow((-1), (- kcolumn/2));
                if (abs(me_E2) < pow(10, -10)) { me_E2 = 0; }
            }
            else if (kcolumn%2 == 0){
                me_E2 = + 2*wigner_d_small_pihalf(jrow)(jrow + kcolumn, jrow + krow)*cos(((kcolumn + krow)*M_PI/4 + 2*M_PI/3))
                        + 2*wigner_d_small_pihalf(jrow)(jrow - kcolumn, jrow + krow)*cos(((- kcolumn + krow)*M_PI/4 + 2*M_PI/3))*pow((-1), jrow)*pow((-1), (- kcolumn/2));
                if (abs(me_E2) < pow(10, -10)) { me_E2 = 0; }
            }
            else {
                me_E2 = 0;
            }
        }
        else {
            me_E2 = 0;
        }
    }
    else {
        me_E2 = 0;
    }
    return me_E2;
}

double SymmetrizeMatrix::F1() {
    double me_F1;
    if (krow%2 == 0) {
        if (jrow == jcolumn && n0row == n0column && nplsrow == nplscolumn && nmnsrow == nmnsrow && mBrow == mBcolumn) {
            if (krow == 0 && kcolumn == 0) {
                me_F1 = 1 - pow((-1), jrow)*pow((-1), (krow/2));
            if (abs(me_F1) < pow(10, -10)) { me_F1 = 0; }
            }
            else if (krow == kcolumn) {
                me_F1 = 1;
            if (abs(me_F1) < pow(10, -10)) { me_F1 = 0; }
            }
            else if (krow == - kcolumn) {
                me_F1 = - pow((-1), jrow)*pow((-1), (krow/2));
            if (abs(me_F1) < pow(10, -10)) { me_F1 = 0; }
            }
            else {
                me_F1 = 0;
            }
        }
        else {
            me_F1 = 0;
        }
    }
    else {
        me_F1 = 0;
    }
    return me_F1;
}

complex<double> SymmetrizeMatrix::F2() {
    complex i(0.0, 1.0);
    complex<double> me_F2;
    if (krow%2 == 1 || krow%2 == -1) {
        if (jrow == jcolumn && n0row == n0column && nplsrow == nplscolumn && nmnsrow == nmnsrow && mBrow == mBcolumn) {
            if (kcolumn%2 == 1 || kcolumn%2 == -1) {
                me_F2 = + pow((- 1), (kcolumn + krow))*wigner_d_small_pihalf(jrow)(jrow + kcolumn, jrow + krow)*exp((kcolumn + krow)*M_PI/4*i)
                        + pow((- 1), (- kcolumn + krow))*wigner_d_small_pihalf(jrow)(jrow - kcolumn, jrow + krow)*exp((- kcolumn + krow)*M_PI/4*i)*i*pow((-1), jrow)*pow((-1), ((- kcolumn - 1)/2));
                if (abs(me_F2.real()) < pow(10, -10)) { me_F2.real(0); }
                if (abs(me_F2.imag()) < pow(10, -10)) { me_F2.imag(0); }
            }
            else {
                me_F2 = 0;
            }
        }
        else {
            me_F2 = 0;
        }
    }
    else {
        me_F2 = 0;
    }
    return me_F2;
}

complex<double> SymmetrizeMatrix::F3() {
    complex i(0.0, 1.0);
    complex<double> me_F3;
    if (krow%2 == 1 || krow%2 == -1) {
        if (jrow == jcolumn && n0row == n0column && nplsrow == nplscolumn && nmnsrow == nmnsrow && mBrow == mBcolumn) {
            if (kcolumn%2 == 1 || kcolumn%2 == -1) {
                me_F3 = + wigner_d_small_pihalf(jrow)(jrow + kcolumn, jrow + krow)*exp(- ( kcolumn + krow)*M_PI/4*(i))
                        + wigner_d_small_pihalf(jrow)(jrow - kcolumn, jrow + krow)*exp(- (- kcolumn + krow)*M_PI/4*(i))*(- i)*pow((-1), jrow)*pow((-1), ((- kcolumn - 1)/2));
                if (abs(me_F3.real()) < pow(10, -10)) { me_F3.real(0); }
                if (abs(me_F3.imag()) < pow(10, -10)) { me_F3.imag(0); }
            }
            else {
                me_F3 = 0;
            }
        }
        else {
            me_F3 = 0;
        }
    }
    else {
        me_F3 = 0;
    }
    return me_F3;
}

vector<vector<complex<double>>> U (int m, string sym) 
{    
    MatrixParameter dim(m, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
    vector<vector<complex<double>>> result(dim.matrix_size(), vector<complex<double>>(dim.matrix_size()));
    vector<int> drow(1, 0);
    vector<int> ddrow(1, 0);
    int dddrow= 0;
    for (int mBrow = - j; mBrow <= j; mBrow++) {                                                    // assign row mB
        for (int jrow = abs(mBrow); jrow <= j; jrow++) {                                            // assign row j
            for (int krow = - jrow; krow <= jrow; krow++) {                                         // assign row k
                for (int n0row = 0; n0row <= n0; n0row++) {                                         // assign row n0
                    for (int vrow = abs(m - mBrow); vrow <= (npls + nmns - abs(m - mBrow)); vrow += 2) {  //assign row v
                        
                        int jrdim = 0;
                        int jdrow;
                        for (int num = 0; num<= j; num++) {
                            jdrow = (2 * abs(num) + 1);
                            jrdim += jdrow;
                        }
                        int mrminus = 0;
                        int mmmrow;
                        for (int mm = 0; mm <= abs(mBrow); mm++) {
                            if (mm == 0) {
                                mmmrow = 0;
                            }
                            else {
                                mmmrow = (2 * (abs(mm) - 1) + 1);
                            }
                            mrminus += mmmrow;
                        }
                        int d = (npls + 1 - abs(m - mBrow))*(jrdim - mrminus);
                        drow.push_back(d);
                        if (d != drow[drow.size() - 2]) {
                            ddrow.push_back(d);
                            dddrow += ddrow[ddrow.size() - 2];
                        }
                        vector<int> dcolumn(1, 0);
                        vector<int> ddcolumn(1, 0);
                        int dddcolumn = 0;
                        for (int mBcolumn = - j; mBcolumn <= j; mBcolumn++) {                                                           // assign column mB
                            for (int jcolumn = abs(mBcolumn); jcolumn <= j; jcolumn++) {                                                // assign column j      
                                for (int kcolumn = - jcolumn; kcolumn <= jcolumn; kcolumn++) {                                          // assign column k
                                    for (int n0column = 0; n0column <= n0; n0column++) {                                                // assign column n0
                                        for (int vcolumn = abs(m - mBcolumn); vcolumn <= (npls + nmns - abs(m - mBcolumn)); vcolumn += 2) {   // assign column v
                                            int jcdim = 0;
                                            int jdcolumn;
                                            for (int num = 0; num<= j; num++) {
                                                jdcolumn = (2 * abs(num) + 1);
                                                jcdim += jdcolumn;
                                                }
                                            int mcminus = 0;
                                            int mmmcolumn;
                                            for (int mm = 0; mm <= abs(mBcolumn); mm++) {
                                                if (mm == 0) {
                                                    mmmcolumn = 0;
                                                }
                                                else {
                                                    mmmcolumn = (2 * (abs(mm) - 1) + 1);
                                                }
                                                mcminus += mmmcolumn;
                                            }
                                            int d = (npls + 1 - abs(m - mBcolumn))*(jcdim - mcminus);
                                            dcolumn.push_back(d);
                                            if (d != dcolumn[dcolumn.size() - 2]) {
                                                ddcolumn.push_back(d);
                                                dddcolumn += ddcolumn[ddcolumn.size() - 2];
                                            }
                                            // matrix element
                                            SymmetrizeMatrix matrix(m, jrow, mBrow, krow, n0row, vrow, jcolumn, mBcolumn, kcolumn, n0column, vcolumn);
                                            if      (sym == "A1") { result[dddrow + matrix.row()][dddcolumn + matrix.column()] = matrix.A1(); }
                                            else if (sym == "E1") { result[dddrow + matrix.row()][dddcolumn + matrix.column()] = matrix.E1(); }
                                            else if (sym == "E2") { result[dddrow + matrix.row()][dddcolumn + matrix.column()] = matrix.E2(); }
                                            else if (sym == "F1") { result[dddrow + matrix.row()][dddcolumn + matrix.column()] = matrix.F1(); }
                                            else if (sym == "F2") { result[dddrow + matrix.row()][dddcolumn + matrix.column()] = matrix.F2(); }
                                            else if (sym == "F3") { result[dddrow + matrix.row()][dddcolumn + matrix.column()] = matrix.F3(); }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    return result;
}

tuple<Eigen::Vector<vector<vector<complex<double>>>, 6>, Eigen::Vector<vector<vector<complex<double>>>, 3>, Eigen::Vector<vector<vector<complex<double>>>, 1>> Projects (int m)
{
    
    Eigen::Vector<vector<vector<complex<double>>>, 6> projects6;

    clock_t start_A1 = clock();
    projects6(0) = U(m, "A1");
    clock_t end_A1 = clock();

    clock_t start_E1 = clock();
    projects6(1) = U(m, "E1");
    clock_t end_E1 = clock();

    clock_t start_E2 = clock();
    projects6(2) = U(m, "E2");
    clock_t end_E2 = clock();

    clock_t start_F1 = clock();
    projects6(3) = U(m, "F1");
    clock_t end_F1 = clock();

    clock_t start_F2 = clock();
    projects6(4) = U(m, "F2");
    clock_t end_F2 = clock();

    clock_t start_F3 = clock();
    projects6(5) = U(m, "F3");
    clock_t end_F3 = clock();

    const double time_A1 = static_cast<double>(end_A1 - start_A1) / CLOCKS_PER_SEC;
    printf("A1 matrix time %lf[s]\n", time_A1);
    const double time_E1 = static_cast<double>(end_E1 - start_E1) / CLOCKS_PER_SEC;
    printf("E1 matrix time %lf[s]\n", time_E1);
    const double time_E2 = static_cast<double>(end_E2 - start_E2) / CLOCKS_PER_SEC;
    printf("E2 matrix time %lf[s]\n", time_E2);
    const double time_F1 = static_cast<double>(end_F1 - start_F1) / CLOCKS_PER_SEC;
    printf("F1 matrix time %lf[s]\n", time_F1);
    const double time_F2 = static_cast<double>(end_F2 - start_F2) / CLOCKS_PER_SEC;
    printf("F2 matrix time %lf[s]\n", time_F2);
    const double time_F3 = static_cast<double>(end_F3 - start_F3) / CLOCKS_PER_SEC;
    printf("F3 matrix time %lf[s]\n", time_F3);

    ofstream time("result_run_time.txt", ios::app);    //txtファイル書き出し
    time << "A1 matrix time " << time_A1 << "[s]\n";
    time << "E1 matrix time " << time_E1 << "[s]\n";
    time << "E2 matrix time " << time_E2 << "[s]\n";
    time << "F1 matrix time " << time_F1 << "[s]\n";
    time << "F2 matrix time " << time_F2 << "[s]\n";
    time << "F3 matrix time " << time_F3 << "[s]\n";
    time.close();

    Eigen::Vector<vector<vector<complex<double>>>, 3> projects3;
    Eigen::Vector<vector<vector<complex<double>>>, 1> projects1;
    

int p = 0;
for (auto& project : projects6) {
    clock_t start_sym = clock();
    p += 1;
    MatrixParameter dim(m, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);

    //remove zero vector
    vector<complex<double>> zero(dim.matrix_size(), 0);
    int rrow = -1;
    for (int row = 0; row < dim.matrix_size(); row++) {
            rrow += 1;
            if (project[rrow] == zero)
			{
			project.erase(project.begin()+rrow);
            rrow -= 1;
			}    
    }

    //normalize eigen vector
    for (int rr = 0; rr < project.size(); rr++) {
        double num = 0;
        double maxindex = 0;
        //normalize eigen vector(each row)
        vector<double> project_abs;
        for (int cc = 0; cc < project[rr].size(); cc++) {
            num += abs(project[rr][cc])*abs(project[rr][cc]);
            project_abs.push_back(abs(project[rr][cc]));
        }
        vector<double>::iterator iter = max_element(project_abs.begin(), project_abs.end());
        size_t index = distance(project_abs.begin(), iter);	

        complex<double> phase = project[rr][index]/abs(project[rr][index]);
        for (int cc = 0; cc < project[rr].size(); cc++) {
            if ( p == 2) {
            }	
            project[rr][cc] = project[rr][cc]/(sqrt(num)*phase);
        }
    }		

    //remove duplicate eigen vector
    vector<int> samerow;
    //collect duplicate vector's row number
    for (int rrr = 0; rrr < project.size(); rrr++) {
        for (int rrrr = rrr; rrrr < project.size(); rrrr++) {
            int bools_sum = 0;
            if (rrr < rrrr) {
                for (int ccc = 0;  ccc < project[rrrr].size(); ccc++) {
                    if (abs(project[rrrr][ccc].real() - project[rrr][ccc].real()) < pow(10, -6) && abs(project[rrrr][ccc].imag() - project[rrr][ccc].imag()) < pow(10, -6)) {
                        bools_sum += 1;
                    }
                    else if (abs(project[rrrr][ccc].real() + project[rrr][ccc].real()) < pow(10, -6) && abs(project[rrrr][ccc].imag() + project[rrr][ccc].imag()) < pow(10, -6)) {
                        bools_sum += 1;
                    }
                }
            }
            if (bools_sum == project[0].size()) {
                samerow.push_back(rrrr);
            }
        }
    }

    sort(samerow.begin(), samerow.end()); 
    samerow.erase(unique(samerow.begin(), samerow.end()), samerow.end()); 

    int samerow_length = samerow.size();
    for (int component = 0; component < samerow_length; component++){
        if (samerow.size() != 0) {
            project.erase(project.begin() + samerow[samerow.size() - 1]);
            samerow.erase(samerow.begin() + (samerow.size() - 1));
        }
    }
    

    if (p == 1) {
        for (int r = 0; r < project.size(); r++) {
            projects3(0).push_back(project[r]);
            projects1(0).push_back(project[r]);
            }
        
    }
    else if (p == 2 || p == 3) {
        for (int r = 0; r < project.size(); r++) {
            projects3(1).push_back(project[r]);
            projects1(0).push_back(project[r]);
            }
    }
    else if (p == 4 || p == 5 || p == 6) {
        for (int r = 0; r < project.size(); r++) {
            projects3(2).push_back(project[r]);
            projects1(0).push_back(project[r]);
            }
    }
    
    clock_t end_sym = clock();
    double time_sym = double(end_sym - start_sym) / CLOCKS_PER_SEC;
    ofstream time2("result_run_time.txt", ios::app);    //txtファイル書き出し
    if      (p == 1){printf("A1 projection time %lf[s]\n", time_sym); time2 << "A1 projection time " << time_sym << "[s]\n";}
    else if (p == 2){printf("E1 projection time %lf[s]\n", time_sym); time2 << "E1 projection time " << time_sym << "[s]\n";}
    else if (p == 3){printf("E2 projection time %lf[s]\n", time_sym); time2 << "E2 projection time " << time_sym << "[s]\n";}
    else if (p == 4){printf("F1 projection time %lf[s]\n", time_sym); time2 << "F1 projection time " << time_sym << "[s]\n";}
    else if (p == 5){printf("F2 projection time %lf[s]\n", time_sym); time2 << "F2 projection time " << time_sym << "[s]\n";}
    else if (p == 6){printf("F3 projection time %lf[s]\n", time_sym); time2 << "F3 projection time " << time_sym << "[s]\n";}
    time2.close();

}

    return make_tuple(projects6, projects3, projects1);
}


Eigen::Vector<vector<vector<complex<double>>>, 4> Projects6 (int m)
{
    
    Eigen::Vector<vector<vector<complex<double>>>, 4> projects6;

    clock_t start_A1 = clock();
    projects6(0) = U(m, "A1");
    clock_t end_A1 = clock();

    clock_t start_E1 = clock();
    projects6(1) = U(m, "E1");
    clock_t end_E1 = clock();

    clock_t start_E2 = clock();
    projects6(2) = U(m, "E2");
    clock_t end_E2 = clock();

    clock_t start_F1 = clock();
    projects6(3) = U(m, "F1");
    clock_t end_F1 = clock();
    

    const double time_A1 = static_cast<double>(end_A1 - start_A1) / CLOCKS_PER_SEC;
    printf("A1 matrix time %lf[s]\n", time_A1);
    const double time_E1 = static_cast<double>(end_E1 - start_E1) / CLOCKS_PER_SEC;
    printf("E1 matrix time %lf[s]\n", time_E1);
    const double time_E2 = static_cast<double>(end_E2 - start_E2) / CLOCKS_PER_SEC;
    printf("E2 matrix time %lf[s]\n", time_E2);
    const double time_F1 = static_cast<double>(end_F1 - start_F1) / CLOCKS_PER_SEC;
    printf("F1 matrix time %lf[s]\n", time_F1);


    ofstream time("result_run_time.txt", ios::app);    //txtファイル書き出し
    time << "A1 matrix time " << time_A1 << "[s]\n";
    time << "E1 matrix time " << time_E1 << "[s]\n";
    time << "E2 matrix time " << time_E2 << "[s]\n";
    time << "F1 matrix time " << time_F1 << "[s]\n";
    time.close();

int p = 0;
for (auto& project : projects6) {
    clock_t start_sym = clock();
    p += 1;
    MatrixParameter dim(m, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);

    //remove zero vector
    vector<complex<double>> zero(dim.matrix_size(), 0);
    int rrow = -1;
    for (int row = 0; row < dim.matrix_size(); row++) {
            rrow += 1;
            if (project[rrow] == zero)
			{
			project.erase(project.begin()+rrow);
            rrow -= 1;
			}    
    }

    //normalize eigen vector
    for (int rr = 0; rr < project.size(); rr++) {
        double num = 0;
        double maxindex = 0;
        //normalize eigen vector(each row)
        vector<double> project_abs;
        for (int cc = 0; cc < project[rr].size(); cc++) {
            num += abs(project[rr][cc])*abs(project[rr][cc]);
            project_abs.push_back(abs(project[rr][cc]));
        }
        vector<double>::iterator iter = max_element(project_abs.begin(), project_abs.end());
        size_t index = distance(project_abs.begin(), iter);	

        complex<double> phase = project[rr][index]/abs(project[rr][index]);
        for (int cc = 0; cc < project[rr].size(); cc++) {
            project[rr][cc] = project[rr][cc]/(sqrt(num)*phase);
        }
    }		

    //remove duplicate eigen vector
    vector<int> samerow;
    //collect duplicate vector's row number
    for (int rrr = 0; rrr < project.size(); rrr++) {
        for (int rrrr = rrr; rrrr < project.size(); rrrr++) {
            int bools_sum = 0;
            if (rrr < rrrr) {
                for (int ccc = 0;  ccc < project[rrrr].size(); ccc++) {
                    if (abs(project[rrrr][ccc].real() - project[rrr][ccc].real()) < pow(10, -6) && abs(project[rrrr][ccc].imag() - project[rrr][ccc].imag()) < pow(10, -6)) {
                        bools_sum += 1;
                    }
                    else if (abs(project[rrrr][ccc].real() + project[rrr][ccc].real()) < pow(10, -6) && abs(project[rrrr][ccc].imag() + project[rrr][ccc].imag()) < pow(10, -6)) {  //規格化フェイズで位相反転して規格化していた時のための保険
                        bools_sum += 1;
                    }
                }
            }
            if (bools_sum == project[0].size()) {
                samerow.push_back(rrrr);
            }
        }
    }

    sort(samerow.begin(), samerow.end());
    samerow.erase(unique(samerow.begin(), samerow.end()), samerow.end()); 

    int samerow_length = samerow.size();
    for (int component = 0; component < samerow_length; component++){
        if (samerow.size() != 0) {
            project.erase(project.begin() + samerow[samerow.size() - 1]);
            samerow.erase(samerow.begin() + (samerow.size() - 1));
        }
    }

    if (p == 3) {
        for (int r = 0; r < project.size(); r++) {
            projects6(1).push_back(project[r]); 
            }
    }

    clock_t end_sym = clock();
    double time_sym = double(end_sym - start_sym) / CLOCKS_PER_SEC;
    ofstream time2("result_run_time.txt", ios::app); 
    if      (p == 1){printf("A1 projection time %lf[s]\n", time_sym); time2 << "A1 projection time " << time_sym << "[s]\n";}
    else if (p == 2){printf("E1 projection time %lf[s]\n", time_sym); time2 << "E1 projection time " << time_sym << "[s]\n";}
    else if (p == 3){printf("E2 projection time %lf[s]\n", time_sym); time2 << "E2 projection time " << time_sym << "[s]\n";}
    else if (p == 4){printf("F1 projection time %lf[s]\n", time_sym); time2 << "F1 projection time " << time_sym << "[s]\n";}
    time2.close();
}

    return projects6;
}
