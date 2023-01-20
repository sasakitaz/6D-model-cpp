#pragma once
#include <iostream>
#include <vector>
#include <cstdlib>  // abs() for integer
#include <numeric>  //accumulate
#include <bits/stdc++.h>    //pi: M_PI
#include "Wigner3jSymbol.h" //wigner 3j symbol
#include "tools.h"  //Wigner3jSymbol product, DVR
#define PRINT_MAT(X) cout << #X << ":\n" << X << endl << endl

using namespace std;

extern Eigen::MatrixXd Ri;          //from DVR(kzz, az) in tools.h
extern Eigen::MatrixXd Ri_vec;      //from DVR(kzz, az) in tools.h

class MatrixElement : public tools
{
protected:
    int p;   //spherical harmonics parameter
    int q;   //spherical harmonics parameter

public:
    //基底クラスのコンストラクタ呼び出し
    MatrixElement(int temm, int temjrow, int temmBrow, int temkrow, int temn0row, int temvrow, int temjcolumn, int temmBcolumn, int temkcolumn, int temn0column, int temvcolumn) : tools(temm, temjrow, temmBrow, temkrow, temn0row, temvrow, temjcolumn, temmBcolumn, temkcolumn, temn0column, temvcolumn)
    {
    }
    double Ts();    //stretch kitnetic
    double Tb();    //bend kinetic
    double Tr();    //internal rotation kinetic
    double Vs(double kzz, double az);   //stretch potential
    double Vb();    //bend potential
    double Vr(int p, int q);    //internal rotation potential 
    double Vsb(double kzz, double az);   //stretch-bend coupling potential
    double Vsr(int p, int q);   //stretch-internal rotation coupling potential
    double Vssr(int p, int q);  //stretch-internal rotation coupling potential(stretch Taylor expansion higher term)
    double Vbr(int p, int q);   //bend-internal rotation coupling potential
    double angular_coupling(int p, int q, int m_coiple);
    double Vam(int p, int q, int m_coiple); //angular momentum coupling potential
    double VamR(int p, int q, int m_couple);    //angular momentum coupling potential(R dependence)
};

double MatrixElement::Ts(){
    double me_Ts;
    if (jrow == jcolumn && krow == kcolumn && nplsrow == nplscolumn && nmnsrow == nmnscolumn){
        switch (n0row - n0column){
            case 0:
                me_Ts = 2*n0column + 1;
                break;
            case + 2:
                me_Ts = - sqrt((n0column + 1)*(n0column + 2));
                break;
            case - 2:
                me_Ts = - sqrt(n0column*(n0column - 1));
                break;
            default:
                me_Ts = 0;
                break;
        }
    }
    else {
        me_Ts = 0;
    }
    return me_Ts;
}

double MatrixElement::Tb(){
    double me_Tb;
    if (jrow ==jcolumn && krow == kcolumn && n0row == n0column) {
        if (nplsrow - nplscolumn == + 1 && nmnsrow - nmnscolumn == + 1){
            me_Tb = - 2*sqrt((nplscolumn + 1)*(nmnscolumn + 1));
        }
        else if (nplsrow - nplscolumn == - 1 && nmnsrow - nmnscolumn == - 1) {
            me_Tb = - 2*sqrt((nplscolumn)*(nmnscolumn));
        }
        else if (nplsrow - nplscolumn == 0 && nmnsrow - nmnscolumn == 0) {
            me_Tb = 2*(nplscolumn + nmnscolumn + 1);
        }
        else {
            me_Tb = 0;
        }
    }
    else {
        me_Tb = 0;
    }
    return me_Tb;
}

double MatrixElement::Tr(){
    double me_Tr;
    if (jrow == jcolumn && krow == kcolumn && n0row == n0column && nplsrow == nplscolumn && nmnsrow == nmnscolumn) {
        me_Tr = jrow*(jrow + 1);
    }
    else {
        me_Tr = 0;
    }
    return me_Tr;
}

double MatrixElement::Vs(double kzz, double az){
    double me_Vs;
    //PRINT_MAT(Ri);
    double VVs = 0;
    for (int i = 0; i <= n0; i++) {
        VVs += (kzz*pow(1 - (exp(-az*Ri(i, 0))), 2))*Ri_vec(n0row, i)*Ri_vec(n0column, i);
    }
    //Vs
    if (jrow == jcolumn && krow == kcolumn && nplsrow == nplscolumn && nmnsrow == nmnscolumn) {
        me_Vs = VVs;
    }
    else {
        me_Vs = 0;
    }
    return me_Vs;
}

double MatrixElement::Vb() {
    double me_Vb;
    if (jrow == jcolumn && krow == kcolumn && n0row == n0column) {
        if (nplsrow - nplscolumn == + 1 && nmnsrow - nmnscolumn == + 1){
            me_Vb = 2*sqrt((nplscolumn + 1)*(nmnscolumn + 1));
        }
        else if (nplsrow - nplscolumn == - 1 && nmnsrow - nmnscolumn == - 1) {
            me_Vb = 2*sqrt((nplscolumn)*(nmnscolumn));
        }
        else if (nplsrow - nplscolumn == 0 && nmnsrow - nmnscolumn == 0) {
            me_Vb = 2*(nplscolumn + nmnscolumn + 1);
        }
        else {
            me_Vb = 0;
        }
    }
    else {
        me_Vb = 0;
    }
    return me_Vb;
}

double MatrixElement::Vr(int p, int q) {
    double me_Vr;
    if (n0row == n0column && nplsrow == nplscolumn && nmnsrow == nmnscolumn) {
        me_Vr = pow((-1), (mBrow - kcolumn))*sqrt((2*jcolumn + 1)*(2*jrow + 1))*tools::Wigner3jSymbol(p, q);
    }
    else {
        me_Vr = 0;
    }
    return me_Vr;
}

double MatrixElement::Vsb(double kzz, double az) {
    //unitary matrix product
    double VVsb = 0;
    for (int i = 0; i <= n0; i++) {
        VVsb += kxxz*(1 - exp(-az*Ri(i, 0)))*Ri_vec(n0row, i)*Ri_vec(n0column, i);
    }
    //cout << VVsb << endl;
    //Vsb
    double me_Vsb;
    if (jrow == jcolumn && krow == kcolumn) {
        if (nplsrow == nplscolumn + 1 && nmnsrow == nmnscolumn + 1) {
            me_Vsb = 2*sqrt((nplscolumn + 1)*(nmnscolumn + 1))*VVsb;
        }
        else if (nplsrow == nplscolumn - 1 && nmnsrow == nmnscolumn - 1) {
            me_Vsb = 2*sqrt((nplscolumn)*(nmnscolumn))*VVsb;
        }
        else if (nplsrow - nplscolumn == 0 && nmnsrow - nmnscolumn == 0) {
            me_Vsb = 2*(nplscolumn + nmnscolumn + 1)*VVsb;
        }
        else {
            me_Vsb = 0;
        }
    }
    else {
        me_Vsb = 0;
    }
    return me_Vsb;
}

double MatrixElement::Vsr(int p, int q) {
    double me_Vsr;
    if (nplsrow == nplscolumn && nmnsrow == nmnscolumn) {
        if (n0row - n0column == + 1) {
            me_Vsr = sqrt(n0column + 1)
                    *pow((-1), (mBcolumn - kcolumn))*sqrt((2*jcolumn + 1)*(2*jrow + 1))*tools::Wigner3jSymbol(p, q);
        }
        else if (n0row - n0column == - 1) {
            me_Vsr = sqrt(n0column)
                    *pow((-1), (mBcolumn - kcolumn))*sqrt((2*jcolumn + 1)*(2*jrow + 1))*tools::Wigner3jSymbol(p, q);
        }
        else {
            me_Vsr = 0;
        }
    }
    else {
        me_Vsr = 0;
    }
    return me_Vsr;
}

double MatrixElement::Vssr(int p, int q) {
    double me_Vssr;
    if (nplsrow == nplscolumn && nmnsrow == nmnscolumn) {
        if (n0row - n0column == 0) {
            me_Vssr = (2*n0column + 1)
                    *pow((-1), (mBcolumn - kcolumn))*sqrt((2*jcolumn + 1)*(2*jrow + 1))*tools::Wigner3jSymbol(p, q);
        }
        else if (n0row - n0column == + 2) {
            me_Vssr = sqrt((n0column + 1)*(n0column + 2))
                    *pow((-1), (mBcolumn - kcolumn))*sqrt((2*jcolumn + 1)*(2*jrow + 1))*tools::Wigner3jSymbol(p, q);
        }
        else if (n0row - n0column == - 2) {
            me_Vssr = sqrt(n0column*(n0column - 1))
                    *pow((-1), (mBcolumn - kcolumn))*sqrt((2*jcolumn + 1)*(2*jrow + 1))*tools::Wigner3jSymbol(p, q);
        }
        else {
            me_Vssr = 0;
        }
    }
    else {
        me_Vssr = 0;
    }
    return me_Vssr;
}

double MatrixElement::Vbr(int p, int q) {
    double me_Vbr;
    if (n0row == n0column) {
        if (nplsrow - nplscolumn == + 1 && nmnsrow - nmnscolumn == + 1) {
            me_Vbr = 2*sqrt((nplscolumn + 1)*(nmnscolumn + 1))
                        *pow((-1), (mBrow - krow))*sqrt((2*jcolumn + 1)*(2*jrow + 1))*tools::Wigner3jSymbol(p, q);
        }
        else if (nplsrow - nplscolumn == - 1 && nmnsrow - nmnscolumn == - 1) {
            me_Vbr = 2*sqrt((nplscolumn)*(nmnscolumn))
                        *pow((-1), (mBrow - krow))*sqrt((2*jcolumn + 1)*(2*jrow + 1))*tools::Wigner3jSymbol(p, q);
        }
        else if (nplsrow - nplscolumn == 0 && nmnsrow - nmnscolumn == 0) {
            me_Vbr = 2*(nplscolumn + nmnscolumn + 1)
                        *pow((-1), (mBrow - krow))*sqrt((2*jcolumn + 1)*(2*jrow + 1))*tools::Wigner3jSymbol(p, q);
        }
        else {
            me_Vbr = 0;
        }
    }
    else {
        me_Vbr = 0;
    }
    return me_Vbr;
}

double MatrixElement::angular_coupling(int p, int q, int m_couple) {
    double me_Vam = 0;
    if (m_couple == + 1 ) {
        if (nplsrow - nplscolumn == + 1 && nmnsrow - nmnscolumn == 0) {
            me_Vam =  sqrt(nplscolumn + 1) 
                      *pow((-1), (mBcolumn - kcolumn))*sqrt((2*jcolumn + 1)*(2*jrow + 1))*tools::Wigner3jSymbol_couple(p, q, - 1);
        }
        else if (nplsrow - nplscolumn == 0 && nmnsrow - nmnscolumn == - 1) {
            me_Vam =  sqrt(nmnscolumn)
                      *pow((-1), (mBcolumn - kcolumn))*sqrt((2*jcolumn + 1)*(2*jrow + 1))*tools::Wigner3jSymbol_couple(p, q, - 1);
        }
        else {
            me_Vam = 0;
        }
    }

    else if (m_couple == - 1) {
        if (nplsrow - nplscolumn == - 1 && nmnsrow - nmnscolumn == 0) {
            me_Vam =  sqrt(nplscolumn) 
                      *pow((-1), (mBcolumn - kcolumn))*sqrt((2*jcolumn + 1)*(2*jrow + 1))*tools::Wigner3jSymbol_couple(p, q, + 1);
        }
        else if (nplsrow - nplscolumn == 0 && nmnsrow - nmnscolumn == + 1) {
            me_Vam =  sqrt(nmnscolumn + 1)
                      *pow((-1), (mBcolumn - kcolumn))*sqrt((2*jcolumn + 1)*(2*jrow + 1))*tools::Wigner3jSymbol_couple(p, q, + 1);
        }
        else {
            me_Vam = 0;
        }
    }
    else {
        cout << "!! ERROR, m_couple is not 1 !!" << endl;
        return 0;
    }
    return me_Vam;
}

double MatrixElement::Vam(int p, int q, int m_couple) {
    double me_Vam = 0;
    if (n0row - n0column == 0) {
        me_Vam = MatrixElement::angular_coupling(p, q, m_couple);
    }
    else {
        me_Vam= 0;
    }
    return me_Vam;
}

double MatrixElement::VamR(int p, int q, int m_couple) {
    double VVamR = 0;
    for (int i = 0; i <= n0; i++) {
        VVamR += pow((Ri(i, 0) + Re), (-n))*Ri_vec(n0row, i)*Ri_vec(n0column, i);
    }
    double me_VamR = 0;
    me_VamR = VVamR*MatrixElement::angular_coupling(p, q, m_couple);
    return me_VamR;
}

Eigen::MatrixXd H_H (int m) 
{    
    MatrixParameter dim(m, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
    Eigen::MatrixXd result = Eigen::MatrixXd::Ones(dim.matrix_size(), dim.matrix_size());
    vector<int> drow(1, 0);
    vector<int> ddrow(1, 0);
    int dddrow= 0;
    for (int mBrow = - j; mBrow <= j; mBrow++) {                                                          // assign row mB
        for (int jrow = abs(mBrow); jrow <= j; jrow++) {                                                  // assign row j
            for (int krow = - jrow; krow <= jrow; krow++) {                                               // assign row k
                for (int n0row = 0; n0row <= n0; n0row++) {                                               // assign row n0
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
                                            MatrixElement matrix(m, jrow, mBrow, krow, n0row, vrow, jcolumn, mBcolumn, kcolumn, n0column, vcolumn);
                                            result(dddrow + matrix.row(), dddcolumn + matrix.column())
                                             = 
                                             
                                              + shbaromega/4*matrix.Ts()    //Ts
                                              + matrix.Vs(kzz, az)  //Vs
                                            
                                              + bhbaromega/4*matrix.Tb()    //Tb
                                              + bhbaromega/4*matrix.Vb()    //Vb
                                            
                                              + B*matrix.Tr()   //Tr

                                              //Vr
                                              + V3*( matrix.Vr(3, 2) + matrix.Vr(3, -2) )
                                              + V4*( sqrt(14)*matrix.Vr(4, 0) - sqrt(5)*(matrix.Vr(4, 4) + matrix.Vr(4, -4)) )
                                            
                                              //Vsb
                                              + bhbaromega/(4*kxx)*matrix.Vsb(kzz, az)
                                              //Vsr
                                              + VR3*sqrt( shbaromega/(4*kzz*az*az) )*( matrix.Vsr(3, 2) + matrix.Vsr(3, -2) )
                                              + VR4*sqrt( shbaromega/(4*kzz*az*az) )*( sqrt(14)*matrix.Vsr(4, 0) - sqrt(5)*(matrix.Vsr(4, 4) + matrix.Vsr(4, -4)) )
                                              //Vssr
                                              + VRR3*( shbaromega/(4*kzz*az*az) )*( matrix.Vssr(3, 2) + matrix.Vssr(3, -2) )
                                              + VRR4*( shbaromega/(4*kzz*az*az) )*( sqrt(14)*matrix.Vssr(4, 0) - sqrt(5)*(matrix.Vssr(4, 4) + matrix.Vssr(4, -4)) )
                                              //Vbr
                                              + Vrho3*bhbaromega/(4*kxx)*( matrix.Vbr(3, 2) + matrix.Vbr(3, -2) )
                                              + Vrho4*bhbaromega/(4*kxx)*( sqrt(14)*matrix.Vbr(4, 0) - sqrt(5)*(matrix.Vbr(4, 4) + matrix.Vbr(4, -4)) )
                                              //Vam
                                              + Vam31*bhbaromega/(4*kxx)*( (matrix.Vam(3, 2, 1) + matrix.Vam(3, 2, -1)) + (matrix.Vam(3, -2, 1) + matrix.Vam(3, -2, -1)) )
                                              //VamR
                                              + VamR31*bhbaromega/(4*kxx)*( (matrix.VamR(3, 2, 1) + matrix.VamR(3, 2, -1)) + (matrix.VamR(3, -2, 1) + matrix.VamR(3, -2, -1)) );
                                              
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