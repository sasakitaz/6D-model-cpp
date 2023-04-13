#pragma once
#include <iostream>
#include <vector>
#include <cstdlib>  // abs() for integer
#include <numeric>  //accumulate
#include <cmath>
#include <bits/stdc++.h>    //pi: M_PI
#include "Eigen/Core"   // x(row, column)
#include "Eigen/Dense"  
#define PRINT_MAT(X) cout << #X << ":\n" << X << endl << endl   //output matrix form

using namespace std;

//double M_PI = 3.141592653589793;

//global parameter
//frequantly changeing parameter
const extern __declspec(selectany) int calcm = 0;  //good quantum number m
const extern __declspec(selectany) double diag_coeff = -1;   //Minimum value of eigenvector coefficients to be written to vector.txt file
// matrix size
const extern __declspec(selectany) int n0 = 1;
const extern __declspec(selectany) int npls = 2;
const extern __declspec(selectany) int nmns = 2;
const extern __declspec(selectany) int j = 3;   

//residuals < 50
//potential parameter
double V3 = -50.113;
double V4 = 2.8980;

double VR3 = -191.57;
double VR4 = 11.360;

double VRR3 = 322.19;
double VRR4 = 0;

double Vrho3 = 42.15;
double Vrho4 = -1.0367;

double Vam31 = 477.34;
double VamR31 = - 2173120;

double kxx = 145.14;
double kzz = 622.44;
double kxxz= -183.70;
double az = 1.4440;
double kxxxx = 0;

double n = 6;
double Re = 3.6316;

//constant
double hbar = 1.05457*pow(10, -34);
double c = 2.998*pow(10, 8);
double Iz = 3.81036*pow(10, -25);
double Ix = 1.45451*pow(10, -25);
double B = 5.241;
double mu = 2.20541*pow(10, -26);
double mub = pow(1/mu + Re*Re/Ix, -1);
double omegas = 1/(2*M_PI)*sqrt(2*kzz*az*az*1.98645*pow(10, -23)*pow(10, 20)/mu);
double omegab = 1/(2*M_PI)*sqrt(2*kxx*1.98645*pow(10, -23)*pow(10, 20)/mub);
//double shbaromega = omegas/c/100;
//double bhbaromega = omegab/c/100;
double shbaromega =  10/(2*M_PI*2.998)*sqrt(6.02*1.98645*2/1.328*kzz*az*az); //Bz-CH4 #pseudo diatom molecule's reduced masss (hbar omega)/4
double bhbaromega =  10/(2*M_PI*2.998)*sqrt(6.02*1.98645*2*10*kxx*(1/13.28 + 0.1497));  //Bz-CH4 #pseudo diatom molecule's reduced masss (hbar omega)/4

class MatrixParameter 
{
public:
    MatrixParameter(int temm, int temjrow, int temmBrow, int temkrow, int temn0row, int temvrow, int temjcolumn, int temmBcolumn, int temkcolumn, int temn0column, int temvcolumn)
        : m(temm)
        , jrow(temjrow)
        , mBrow(temmBrow)
        , krow(temkrow)
        , n0row(temn0row)
        , vrow(temvrow)
        , lrow(m - mBrow)
        , jcolumn(temjcolumn)
        , mBcolumn(temmBcolumn)
        , kcolumn(temkcolumn)
        , n0column(temn0column)
        , vcolumn(temvcolumn)
        , lcolumn(m - mBcolumn)
    {
    }
    int matrix_size(); //matrix size
    tuple<vector<int>, vector<int>, vector<int>>size_vec(); 
    int row2();
    int row3();
    int row4();
    int column2();
    int column3();
    int column4();
    int row();
    int column();

protected:
    int m;
    int jrow;
    int mBrow;
    int krow;
    int n0row;
    int vrow;
    int lrow;
    int jcolumn;
    int mBcolumn;
    int kcolumn;
    int n0column;
    int vcolumn;
    int lcolumn;
    int nplsrow = (vrow + lrow)/2;
    int nmnsrow = (vrow - lrow)/2;
    int nplscolumn = (vcolumn + lcolumn)/2;
    int nmnscolumn = (vcolumn - lcolumn)/2;

};

int MatrixParameter::matrix_size()
{
    vector<int> dim(1, 0);
    vector<int> allowl(1, 0);
    vector<int> allowmB(1, 0);

    for (int ll = - nmns; ll <= npls; ll++){
        for (int mB = - j; mB <= j; mB++){

            int jdim = 0;
            for (int num = 0; num <= j; num++){
                int jd = (2*abs(num) + 1);
                jdim += jd;
            }

            int mminus = 0;
            int mmm;
            for (int mm = 0; mm <= abs(mB); mm++){
                if (mm == 0){
                    mmm = 0;
                }
                else {
                    mmm = (2 * (abs(mm) - 1) + 1);
                }
                mminus += mmm;
            }

            if (ll + mB == m){
                int d = (npls + 1 - abs(ll))*(jdim - mminus);
                dim.push_back(d);
                allowl.push_back(ll);
                allowmB.push_back(mB);
            }
        }
    }
    int dimension = accumulate(dim.begin(), dim.end(), 0)*(n0 + 1);
    return dimension;
}

tuple<vector<int>, vector<int>, vector<int>> MatrixParameter::size_vec()
{
    vector<int> dim(1, 0);
    vector<int> allowl(1, 0);
    vector<int> allowmB(1, 0);
    dim.pop_back();
    allowl.pop_back();
    allowmB.pop_back();
    for (int ll = - nmns; ll <= npls; ll++){
        for (int mB = - j; mB <= j; mB++){
            int jdim = 0;
            for (int num = 0; num <= j; num++){
                int jd = (2*abs(num) + 1);
                jdim += jd;
            }

            int mminus = 0;
            int mmm;
            for (int mm = 0; mm <= abs(mB); mm++){
                if (mm == 0){
                    mmm = 0;
                }
                else {
                    mmm = (2*(abs(mm) - 1) + 1);
                }
                mminus += mmm;
            }

            if (ll + mB == m){
                int d = (npls + 1 - abs(ll))*(jdim - mminus);
                dim.push_back(d);
                allowl.push_back(ll);
                allowmB.push_back(mB);
            }
        }
    }

    return forward_as_tuple(dim, allowl, allowmB);
}

int MatrixParameter::row2()
{
    int jrsum = 0;
    int jr;
    for (int numjrow = abs(mBrow); numjrow <= jrow; numjrow++){
        if (numjrow == abs(mBrow)){
            jr = 0;
        }
        else {
            jr = (2*(numjrow - 1) + 1);
        }
        jrsum += jr;
    }
    return jrsum;
}

int MatrixParameter::row3()
{
    int krsum = 0;
    for (int numkrow = - jrow; numkrow <= krow; numkrow++){
        krsum += 1;
    }
    return krsum;
}

int MatrixParameter::row4()
{
    int jjrsum = 0;
    for (int jjrow = abs(mBrow); jjrow <= j; jjrow++){
        jjrsum += 2*jjrow + 1;
    }
    return jjrsum;
}

int MatrixParameter::column2()
{
    int jcsum = 0;
    int jc;
    for (int numjcolumn = abs(mBcolumn); numjcolumn <= jcolumn; numjcolumn++){
        if (numjcolumn == abs(mBcolumn)){
            jc = 0;
        }
        else {
            jc = (2*(numjcolumn - 1) + 1);
        }
        jcsum += jc;
    }
    return jcsum;
}

int MatrixParameter::column3()
{
    int kcsum = 0;
    for (int numkcolumn = - jcolumn; numkcolumn <= kcolumn; numkcolumn++){
        kcsum += 1;
    }
    return kcsum;
}

int MatrixParameter::column4()
{
    int jjcsum = 0;
    for (int jjcolumn = abs(mBcolumn); jjcolumn <= j; jjcolumn++){
        jjcsum += 2*jjcolumn + 1;
    }
    return jjcsum;
}

int MatrixParameter::row(){
    return row2() + (row3() - 1) + (vrow - abs(lrow))/2*row4() + n0row*matrix_size()/(n0 + 1);
}

int MatrixParameter::column(){
    return column2() + (column3() - 1) + (vcolumn - abs(lcolumn))/2*column4() + n0column*matrix_size()/(n0 + 1);
}