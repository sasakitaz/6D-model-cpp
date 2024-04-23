/**************************************************************************
ver.200での変更点
角運動量を考慮して運動エネルギー項を修正
参照：W. Kim et al., J. Chem. Phys. 110, 8461 (1999)

20230703
T_5_2：交換関係を計算しなおし
[A_+^+ A_-^+ - A_; A_-] -> [A_+^+ A_-^+ - A_; A_- + 1]

20230707 ver.300での変更点
Bz-Arを再現するように運動エネルギー項を修正

実験値とはよくあっているが、Vamのsqrt(3/(8*M_PI))がない。
***************************************************************************/

#pragma once
#include <iostream>
#include <vector>
#include <cstdlib>  // abs() for integer
#include <numeric>  //accumulate
#include <bits/stdc++.h>    //pi: M_PI
#include "Wigner3jSymbol.h" //wigner 3j symbol
#include "tools.h"  //Wigner3jSymbol product, DVR, classの継承
#define PRINT_MAT(X) cout << #X << ":\n" << X << endl << endl   //行列の出力

using namespace std;

extern Eigen::MatrixXd Ri;          //from DVR(kzz, az) in tools.h
extern Eigen::MatrixXd Ri_vec;      //from DVR(kzz, az) in tools.h

double ReIxsqrt2mubomegas = Re*pow(10, -10)/Ix*sqrt(2*mub*hbar/omegas);
double ReIxsqrt2muomegas = Re*pow(10, -10)/Ix*sqrt(2*mu*hbar/omegas);
double hbar2Ix = hbar*hbar/Ix*5.03412*pow(10, 22); //1 J = 5.03412 \times 10^22 cm-1
double hbar2Iz = hbar*hbar/Iz*5.03412*pow(10, 22);
double Resqrt2mubhbaromegas = hbar2Ix*Re*pow(10, -10)*sqrt(2*mu*omegas/hbar);
double Resqrt2muhbaromegas = hbar2Ix*Re*pow(10, -10)*sqrt(2*mu*omegas/hbar);
double omegasomegab = shbaromega/bhbaromega;
double Vbcoeff = hbar/(2*(mub*omegab))*pow(10, 20);
double Vscoeff = hbar/(2*mu*omegas)*pow(10, 20);
//double Vscoeff = shbaromega/(4*kzz*az*az);
//double Vbcoeff = bhbaromega*sqrt(1 + Re2Ix*mu)/(4*kxx);

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
    double Tb();    //bend kinetic
    double Tb_1();   //
    double Tb_2();   //
    double Ts();    //stretch kitnetic
    double Ts_1();   //
    double T_3_1();   //
    double T_4_1();   //
    double T_4_2();   //
    double T_5_1();   //
    double T_5_2();   //
    double T_6_1();   //
    double T_6_2();   //
    double T_7_1();   //
    double Tr();    //internal rotation kinetic
    double Tr_1();  //
    double Vs(double kzz, double az);   //stretch potential
    double Vb();    //bend potential
    double Vr(int p, int q);    //internal rotation potential 
    double Vsb(double kxxz, double az);   //stretch-bend coupling potential
    double Vsr(int p, int q);   //stretch-internal rotation coupling potential
    double Vssr(int p, int q);  //stretch-internal rotation coupling potential(stretch Taylor expansion higher term)
    double Vbr(int p, int q);   //bend-internal rotation coupling potential
    double angular_coupling(int p, int q, int m_coiple);
    double Vam(int p, int q, int m_coiple); //angular momentum coupling potential
    double VamR(int p, int q, int m_couple);    //angular momentum coupling potential(R dependence)
};

double MatrixElement::Tb(){
  double me_Tb;
  if (jrow ==jcolumn && krow == kcolumn && n0row == n0column) {
      if (nplsrow - nplscolumn == + 1 && nmnsrow - nmnscolumn == + 1){
          me_Tb = 2*sqrt((nplscolumn + 1)*(nmnscolumn + 1));
      }
      else if (nplsrow - nplscolumn == - 1 && nmnsrow - nmnscolumn == - 1) {
          me_Tb = 2*sqrt((nplscolumn)*(nmnscolumn));
      }
      else if (nplsrow - nplscolumn == 0 && nmnsrow - nmnscolumn == 0) {
          me_Tb = - 2*(nplscolumn + nmnscolumn + 1);
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

double MatrixElement::Tb_1(){
  double me_Tb_1;
  if (jrow ==jcolumn && krow == kcolumn) {
    if (n0row - n0column == + 1) {
      if (nplsrow - nplscolumn == + 1 && nmnsrow - nmnscolumn == + 1){
          me_Tb_1 = 2*sqrt((nplscolumn + 1)*(nmnscolumn + 1))*sqrt(n0column + 1);
      }
      else if (nplsrow - nplscolumn == - 1 && nmnsrow - nmnscolumn == - 1) {
          me_Tb_1 = 2*sqrt((nplscolumn)*(nmnscolumn))*sqrt(n0column + 1);
      }
      else if (nplsrow - nplscolumn == 0 && nmnsrow - nmnscolumn == 0) {
          me_Tb_1 = - 2*(nplscolumn + nmnscolumn + 1)*sqrt(n0column + 1);
      }
      else {
          me_Tb_1 = 0;
      }
    }
    else if (n0row - n0column == - 1) {
      if (nplsrow - nplscolumn == + 1 && nmnsrow - nmnscolumn == + 1){
          me_Tb_1 = 2*sqrt((nplscolumn + 1)*(nmnscolumn + 1))*sqrt(n0column);
      }
      else if (nplsrow - nplscolumn == - 1 && nmnsrow - nmnscolumn == - 1) {
          me_Tb_1 = 2*sqrt((nplscolumn)*(nmnscolumn))*sqrt(n0column);
      }
      else if (nplsrow - nplscolumn == 0 && nmnsrow - nmnscolumn == 0) {
          me_Tb_1 = - 2*(nplscolumn + nmnscolumn + 1)*sqrt(n0column);
      }
      else {
          me_Tb_1 = 0;
      }
    }
    else {
      me_Tb_1 = 0;
    }
  }
  else {
    me_Tb_1 = 0;
  }
  return me_Tb_1;
}

double MatrixElement::Tb_2(){
  double me_Tb_2;
  if (jrow ==jcolumn && krow == kcolumn) {
    if (n0row - n0column == + 2) {
      if (nplsrow - nplscolumn == + 1 && nmnsrow - nmnscolumn == + 1){
          me_Tb_2 = 2*sqrt((nplscolumn + 1)*(nmnscolumn + 1))*sqrt((n0column + 1)*(n0column + 2));
      }
      else if (nplsrow - nplscolumn == - 1 && nmnsrow - nmnscolumn == - 1) {
          me_Tb_2 = 2*sqrt((nplscolumn)*(nmnscolumn))*sqrt((n0column + 1)*(n0column + 2));
      }
      else if (nplsrow - nplscolumn == 0 && nmnsrow - nmnscolumn == 0) {
          me_Tb_2 = - 2*(nplscolumn + nmnscolumn + 1)*sqrt((n0column + 1)*(n0column + 2));
      }
      else {
          me_Tb_2 = 0;
      }
    }
    else if (n0row - n0column == - 2) {
      if (nplsrow - nplscolumn == + 1 && nmnsrow - nmnscolumn == + 1){
          me_Tb_2 = 2*sqrt((nplscolumn + 1)*(nmnscolumn + 1))*sqrt((n0column)*(n0column - 1));
      }
      else if (nplsrow - nplscolumn == - 1 && nmnsrow - nmnscolumn == - 1) {
          me_Tb_2 = 2*sqrt((nplscolumn)*(nmnscolumn))*sqrt((n0column)*(n0column - 1));
      }
      else if (nplsrow - nplscolumn == 0 && nmnsrow - nmnscolumn == 0) {
          me_Tb_2 = - 2*(nplscolumn + nmnscolumn + 1)*sqrt((n0column)*(n0column - 1));
      }
      else {
          me_Tb_2 = 0;
      }
    }
    else if (n0row - n0column == 0) {
      if (nplsrow - nplscolumn == + 1 && nmnsrow - nmnscolumn == + 1){
          me_Tb_2 = 2*sqrt((nplscolumn + 1)*(nmnscolumn + 1))*(2*n0column + 1);
      }
      else if (nplsrow - nplscolumn == - 1 && nmnsrow - nmnscolumn == - 1) {
          me_Tb_2 = 2*sqrt((nplscolumn)*(nmnscolumn))*(2*n0column + 1);
      }
      else if (nplsrow - nplscolumn == 0 && nmnsrow - nmnscolumn == 0) {
          me_Tb_2 = - 2*(nplscolumn + nmnscolumn + 1)*(2*n0column + 1);
      }
      else {
          me_Tb_2 = 0;
      }
    }
    else {
      me_Tb_2 = 0;
    }
  }
  else {
    me_Tb_2 = 0;
  }
  return me_Tb_2;
}

double MatrixElement::Ts(){
    double me_Ts;
    if (jrow == jcolumn && krow == kcolumn && nplsrow == nplscolumn && nmnsrow == nmnscolumn){
        switch (n0row - n0column){
            case 0:
                me_Ts = - 2*n0column - 1;
                break;
            case + 2:
                me_Ts = sqrt((n0column + 1)*(n0column + 2));
                break;
            case - 2:
                me_Ts = sqrt(n0column*(n0column - 1));
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

double MatrixElement::Ts_1(){
  double me_Ts_1;
  if (jrow == jcolumn && krow == kcolumn) {
    if (nplsrow - nplscolumn == + 1 && nmnsrow - nmnscolumn == + 1){
      if (n0row - n0column == + 2) {
        me_Ts_1 = sqrt((n0column + 1)*(n0column + 2))*sqrt((nplscolumn + 1)*(nmnscolumn + 1));
      }
      else if (n0row - n0column == - 2) {
        me_Ts_1 = sqrt(n0column*(n0column - 1))*sqrt((nplscolumn + 1)*(nmnscolumn + 1));
      }
      else if (n0row - n0column == 0) {
        me_Ts_1 = (- 2*n0column - 1)*sqrt((nplscolumn + 1)*(nmnscolumn + 1));
      }
      else {
        me_Ts_1 = 0;
      }
    }
    else if (nplsrow - nplscolumn == - 1 && nmnsrow - nmnscolumn == - 1) {
      if (n0row - n0column == + 2) {
        me_Ts_1 = sqrt((n0column + 1)*(n0column + 2))*sqrt((nplscolumn)*(nmnscolumn));
      }
      else if (n0row - n0column == - 2) {
        me_Ts_1 = sqrt(n0column*(n0column - 1))*sqrt((nplscolumn)*(nmnscolumn));
      }
      else if (n0row - n0column == 0) {
        me_Ts_1 = (- 2*n0column - 1)*sqrt((nplscolumn)*(nmnscolumn));
      }
      else {
        me_Ts_1 = 0;
      }
    }
    else if (nplsrow - nplscolumn == 0 && nmnsrow - nmnscolumn == 0) {
      if (n0row - n0column == + 2) {
        me_Ts_1 = sqrt((n0column + 1)*(n0column + 2))*(nplscolumn + nmnscolumn + 1);
      }
      else if (n0row - n0column == - 2) {
        me_Ts_1 = sqrt(n0column*(n0column - 1))*(nplscolumn + nmnscolumn + 1);
      }
      else if (n0row - n0column == 0) {
        me_Ts_1 = (- 2*n0column - 1)*(nplscolumn + nmnscolumn + 1);
      }
      else {
        me_Ts_1 = 0;
      }
    }
    else {
      me_Ts_1 = 0;
    }
  }
  else {
    me_Ts_1 = 0;
  }
  return me_Ts_1;
}

double MatrixElement::T_3_1(){
  double me_T_3_1;
  if (n0row == n0column && nplsrow == nplscolumn && nmnsrow == nmnscolumn && jrow ==jcolumn && krow == kcolumn) {
    //me_T_3_1 = (4*nplscolumn*(nplscolumn - 1) + 4*nmnscolumn*(nmnscolumn - 1) + 4*nplscolumn*nmnscolumn + 10*nplscolumn + 10*nmnscolumn + 3);
    me_T_3_1 = nplscolumn*nplscolumn + nmnscolumn*nmnscolumn - 2*nplscolumn*nmnscolumn;
  }
  else {
    me_T_3_1 = 0;
  }
  return me_T_3_1;
}

double MatrixElement::T_4_1(){
  double me_T_4_1;
  if (nplsrow == nplscolumn && nmnsrow == nmnscolumn && jrow ==jcolumn && krow == kcolumn) {
    if (n0row - n0column == + 2) {
      me_T_4_1 = sqrt((n0column + 1)*(n0column + 2));
    }
    else if (n0row - n0column == - 2) {
      me_T_4_1 = - sqrt(n0column*(n0column - 1));
    }
    else if (n0row - n0column == 0) {
      me_T_4_1 = 1;
    }
    else {
      me_T_4_1 = 0;
    }
  }
  else {
    me_T_4_1 = 0;
  }
  return me_T_4_1;
}

double MatrixElement::T_5_1(){
  double me_T_5_1;
  if (jrow ==jcolumn && krow == kcolumn) {
    if (n0row - n0column == + 2) {
      if (nplsrow - nplscolumn == + 1 && nmnsrow - nmnscolumn == + 1) {
        me_T_5_1 = sqrt((nplscolumn + 1)*(nmnscolumn + 1))*sqrt((n0column + 1)*(n0column + 2));
      }
      else if (nplsrow - nplscolumn == - 1 && nmnsrow - nmnscolumn == - 1) {
        me_T_5_1 = - sqrt(nplscolumn*nmnscolumn)*sqrt((n0column + 1)*(n0column + 2));
      }
      else if (nplsrow - nplscolumn == 0 && nmnsrow - nmnscolumn == 0) {
        me_T_5_1 = sqrt((n0column + 1)*(n0column + 2));
      }
      else {
        me_T_5_1 = 0;
      }
    }
    else if (n0row - n0column == - 2) {
      if (nplsrow - nplscolumn == + 1 && nmnsrow - nmnscolumn == + 1) {
        me_T_5_1 = - sqrt((nplscolumn + 1)*(nmnscolumn + 1))*sqrt(n0column*(n0column - 1));
      }
      else if (nplsrow - nplscolumn == - 1 && nmnsrow - nmnscolumn == - 1) {
        me_T_5_1 = sqrt(nplscolumn*nmnscolumn)*sqrt(n0column*(n0column - 1));
      }
      else if (nplsrow - nplscolumn == 0 && nmnsrow - nmnscolumn == 0) {
        me_T_5_1 = - sqrt((n0column)*(n0column));
      }
      else {
        me_T_5_1 = 0;
      }
    }
    else {
      me_T_5_1 = 0;
    }
  }
  else {
    me_T_5_1 = 0;
  }
  return me_T_5_1;
}

double MatrixElement::T_5_2(){
  double me_T_5_2;
  if (jrow - jcolumn == 0 && krow - kcolumn == 0) {
    if (n0row - n0column == + 1) {
      if (nplsrow - nplscolumn == + 1 && nmnsrow - nmnscolumn == + 1) {
        me_T_5_2 = sqrt((nplscolumn + 1)*(nmnscolumn + 1))*sqrt(n0column + 1);
      }
      else if (nplsrow - nplscolumn == - 1 && nmnsrow - nmnscolumn == - 1) {
        me_T_5_2 = - sqrt(nplscolumn*nmnscolumn)*sqrt(n0column + 1);
      }
      //else if (nplsrow - nplscolumn == 0 && nmnsrow - nmnscolumn == 0) {
      //  me_T_5_2 = 2*sqrt(n0column + 1);
      //}
      else {
        me_T_5_2 = 0;
      }
    }
    else if (n0row - n0column == - 1) {
      if (nplsrow - nplscolumn == + 1 && nmnsrow - nmnscolumn == + 1) {
        me_T_5_2 = - sqrt((nplscolumn + 1)*(nmnscolumn + 1))*sqrt(n0column);
      }
      else if (nplsrow - nplscolumn == - 1 && nmnsrow - nmnscolumn == - 1) {
        me_T_5_2 = sqrt(nplscolumn*nmnscolumn)*sqrt(n0column);
      }
      //else if (nplsrow - nplscolumn == 0 && nmnsrow - nmnscolumn == 0) {
      //  me_T_5_2 = 2*sqrt(n0column);
      //}
      else {
        me_T_5_2 = 0;
      }
    }
    else {
      me_T_5_2 = 0;
    }
  }
  else {
    me_T_5_2 = 0;
  }
  return me_T_5_2;
}

double MatrixElement::T_6_1(){
  double me_T_6_1;
  if (jrow == jcolumn && krow == kcolumn) {
    if (n0row - n0column == + 1) {
      if (mBrow - mBcolumn == + 1) {
        if (nplsrow - nplscolumn == - 1) {
          me_T_6_1 = - sqrt(nplscolumn)*sqrt((jcolumn - mBcolumn)*(jcolumn + mBcolumn + 1))*sqrt(n0column + 1);
        }
        else if (nmnsrow - nmnscolumn == + 1) {
          me_T_6_1 = sqrt(nmnscolumn + 1)*sqrt((jcolumn - mBcolumn)*(jcolumn + mBcolumn + 1))*sqrt(n0column + 1);
        }
        else {
          me_T_6_1 = 0;
        }
      }
      else if (mBrow - mBcolumn == - 1) {
        if (nplsrow - nplscolumn == + 1) {
          me_T_6_1 = - sqrt(nplscolumn + 1)*sqrt((jcolumn + mBcolumn)*(jcolumn - mBcolumn + 1))*sqrt(n0column + 1);
        }
        else if (nmnsrow - nmnscolumn == - 1) {
          me_T_6_1 = sqrt(nmnscolumn)*sqrt((jcolumn + mBcolumn)*(jcolumn - mBcolumn + 1))*sqrt(n0column + 1);
        }
        else {
          me_T_6_1 = 0;
        }
      }
      else {
        me_T_6_1 = 0;
      }
    }
    else if (n0row - n0column == - 1) {
      if (mBrow - mBcolumn == + 1) {
        if (nplsrow - nplscolumn == - 1) {
          me_T_6_1 = - sqrt(nplscolumn)*sqrt((jcolumn - mBcolumn)*(jcolumn + mBcolumn + 1))*sqrt(n0column);
        }
        else if (nmnsrow - nmnscolumn == + 1) {
          me_T_6_1 = sqrt(nmnscolumn + 1)*sqrt((jcolumn - mBcolumn)*(jcolumn + mBcolumn + 1))*sqrt(n0column);
        }
        else {
          me_T_6_1 = 0;
        }
      }
      else if (mBrow - mBcolumn == - 1) {
        if (nplsrow - nplscolumn == + 1) {
          me_T_6_1 = - sqrt(nplscolumn + 1)*sqrt((jcolumn + mBcolumn)*(jcolumn - mBcolumn + 1))*sqrt(n0column);
        }
        else if (nmnsrow - nmnscolumn == - 1) {
          me_T_6_1 = sqrt(nmnscolumn)*sqrt((jcolumn + mBcolumn)*(jcolumn - mBcolumn + 1))*sqrt(n0column);
        }
        else {
          me_T_6_1 = 0;
        }
      }
      else {
        me_T_6_1 = 0;
      }
    }
    else {
      me_T_6_1 = 0;
    }
  }
  else {
    me_T_6_1 = 0;
  }
  return me_T_6_1;
}

double MatrixElement::T_6_2(){
  double me_T_6_2;
  if (jrow == jcolumn && krow == kcolumn) {
    if (n0row - n0column == + 1) {
      if (mBrow - mBcolumn == + 1) {
        if (nplsrow - nplscolumn == - 1) {
          me_T_6_2 = - sqrt(nplscolumn)*sqrt((jcolumn - mBcolumn)*(jcolumn + mBcolumn + 1))*sqrt(n0column + 1);
        }
        else if (nmnsrow - nmnscolumn == + 1) {
          me_T_6_2 = - sqrt(nmnscolumn + 1)*sqrt((jcolumn - mBcolumn)*(jcolumn + mBcolumn + 1))*sqrt(n0column + 1);
        }
        else {
          me_T_6_2 = 0;
        }
      }
      else if (mBrow - mBcolumn == - 1) {
        if (nplsrow - nplscolumn == + 1) {
          me_T_6_2 = sqrt(nplscolumn + 1)*sqrt((jcolumn + mBcolumn)*(jcolumn - mBcolumn + 1))*sqrt(n0column + 1);
        }
        else if (nmnsrow - nmnscolumn == - 1) {
          me_T_6_2 = sqrt(nmnscolumn)*sqrt((jcolumn + mBcolumn)*(jcolumn - mBcolumn + 1))*sqrt(n0column + 1);
        }
        else {
          me_T_6_2 = 0;
        }
      }
      else {
        me_T_6_2 = 0;
      }
    }
    else if (n0row - n0column == - 1) {
      if (mBrow - mBcolumn == + 1) {
        if (nplsrow - nplscolumn == - 1) {
          me_T_6_2 = sqrt(nplscolumn)*sqrt((jcolumn - mBcolumn)*(jcolumn + mBcolumn + 1))*sqrt(n0column);
        }
        else if (nmnsrow - nmnscolumn == + 1) {
          me_T_6_2 = sqrt(nmnscolumn + 1)*sqrt((jcolumn - mBcolumn)*(jcolumn + mBcolumn + 1))*sqrt(n0column);
        }
        else {
          me_T_6_2 = 0;
        }
      }
      else if (mBrow - mBcolumn == - 1) {
        if (nplsrow - nplscolumn == + 1) {
          me_T_6_2 = - sqrt(nplscolumn + 1)*sqrt((jcolumn + mBcolumn)*(jcolumn - mBcolumn + 1))*sqrt(n0column);
        }
        else if (nmnsrow - nmnscolumn == - 1) {
          me_T_6_2 = - sqrt(nmnscolumn)*sqrt((jcolumn + mBcolumn)*(jcolumn - mBcolumn + 1))*sqrt(n0column);
        }
        else {
          me_T_6_2 = 0;
        }
      }
      else {
        me_T_6_2 = 0;
      }
    }
    else {
      me_T_6_2 = 0;
    }
  }
  else {
    me_T_6_2 = 0;
  }
  return me_T_6_2;
}

double MatrixElement::T_7_1(){
  double me_T_7_1;
  if (jrow == jcolumn && krow == kcolumn && n0row == n0column && nplsrow == nplscolumn && nmnsrow == nmnscolumn) {
    me_T_7_1 = (nplscolumn + nmnscolumn + 1)*mBcolumn;
  }
  else {
    me_T_7_1 = 0;
  }
  return me_T_7_1;
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

double MatrixElement::Tr_1(){
    double me_Tr_1;
    if (jrow == jcolumn && krow == kcolumn && n0row == n0column && nplsrow == nplscolumn && nmnsrow == nmnscolumn) {
        me_Tr_1 = mBcolumn^2;
    }
    else {
        me_Tr_1 = 0;
    }
    return me_Tr_1;
}

double MatrixElement::Vs(double kzz, double az){
    double me_Vs;
    //PRINT_MAT(Ri);
    double VVs = 0;
    for (int i = 0; i <= n0; i++) {
        VVs += (1 - exp(-az*Ri(i, 0)))*(1 - exp(-az*Ri(i, 0)))*Ri_vec(n0row, i)*Ri_vec(n0column, i);
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

double MatrixElement::Vsb(double kxxz, double az) {
    //unitary matrix product
    double VVsb = 0;
    for (int i = 0; i <= n0; i++) {
        VVsb += (1 - exp(-az*Ri(i, 0)))*Ri_vec(n0row, i)*Ri_vec(n0column, i);
    }
    //cout << VVsb << endl;
    //Vsb
    double me_Vsb;
    if (jrow == jcolumn && krow == kcolumn) {
        if (nplsrow - nplscolumn == + 1 && nmnsrow - nmnscolumn == + 1) {
            me_Vsb = 2*sqrt((nplscolumn + 1)*(nmnscolumn + 1))*VVsb;
        }
        else if (nplsrow - nplscolumn == - 1 && nmnsrow - nmnscolumn == - 1) {
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
        /*
        if (nplsrow - nplscolumn == + 1 && nmnsrow - nmnscolumn == + 1) {
            me_Vbr = 2*sqrt((nplscolumn + 1)*(nmnscolumn + 1))
                        *pow((-1), (mBcolumn - kcolumn))*sqrt((2*jcolumn + 1)*(2*jrow + 1))*tools::Wigner3jSymbol(p, q);
        }
        else if (nplsrow - nplscolumn == - 1 && nmnsrow - nmnscolumn == - 1) {
            me_Vbr = 2*sqrt((nplscolumn)*(nmnscolumn))
                        *pow((-1), (mBcolumn - kcolumn))*sqrt((2*jcolumn + 1)*(2*jrow + 1))*tools::Wigner3jSymbol(p, q);
        }
        else if (nplsrow - nplscolumn == 0 && nmnsrow - nmnscolumn == 0) {
            me_Vbr = 2*(nplscolumn + nmnscolumn + 1)
                        *pow((-1), (mBcolumn - kcolumn))*sqrt((2*jcolumn + 1)*(2*jrow + 1))*tools::Wigner3jSymbol(p, q);
        }
        */
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

//Pythonのときの旧プログラムから関数自体を変更
double MatrixElement::angular_coupling(int p, int q, int m_couple) {
    double me_Vam = 0;
    //m_couple = \pm 1のみを想定。それ以外は行列要素が異なるためこの行列要素は使えない。
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
    //if (Wigner3jSymbol_couple(p, q, m_couple) != 0){
    //    cout<< "coodinate, " << vrow - vcolumn << lrow - lcolumn << "         3jsymbol, " << tools::Wigner3jSymbol_couple(p, q, m_couple) << endl;
    //}
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
    cout << "bhbaromega ,\t" << bhbaromega << endl;
    cout << "bhbaromega*sqrt(1 + Re2Ix*mu) ,\t" << bhbaromega*sqrt(1 + Re2Ix*mu) << endl;
    cout << "shbaromega ,\t" << shbaromega << endl;
    cout << "(4*kxx) ,\t" << (4*kxx) << endl;
    cout << "(4*kzz*az*az) ,\t" << (4*kzz*az*az) << endl;

    cout << "hbar2Iz,\t" << hbar*hbar/Iz*5.03412*pow(10, 22) << endl;
    cout << "hbar2Ix,\t" << hbar*hbar/Ix*5.03412*pow(10, 22) << endl;

    cout << "hbar/(mu*omegab) ,\t" << hbar/(mu*omegab) << endl;

    cout << "ReIxsqrt2mubomegas ,\t" << ReIxsqrt2mubomegas << endl;
    cout << "ReIxsqrt2muomegas ,\t" << ReIxsqrt2muomegas << endl;
    cout << "hbar2Ix ,\t" << hbar2Ix << endl;
    cout << "sqrt(1 + Re2Ix*mu) ,\t" << sqrt(1 + Re2Ix*mu) << endl;

    cout << "- Resqrt2mubhbaromegas/2 ,\t" << - Resqrt2mubhbaromegas/2 << endl;

    cout << "matrix.Tb() " << "\t" << bhbaromega/4*(1 + Re2Ix*mu) << endl;
    cout << "matrix.Tb_1()" << "\t" << bhbaromega/4*ReIxsqrt2mubomegas*sqrt(mub/mu) << endl;
    cout << "matrix.Tb_2()" << "\t" << bhbaromega/4*hbar2Ix/2*1/shbaromega*mub/mu << endl;

    cout << "matrix.Ts()" << "\t" << shbaromega/4 << endl;
    cout << "matrix.Ts_1()" << "\t" << hbar2Ix/4*omegas/omegab << "\n" << endl;

    cout << "matrix.T_3_1()" << "\t" << - hbar2Iz << endl;
    cout << "matrix.T_4_1()" << "\t" << + hbar2Ix/2 << endl;
    cout << "matrix.T_5_1()" << "\t" << - hbar2Ix/2 << endl;
    cout << "matrix.T_5_2()" << "\t" << - shbaromega*ReIxsqrt2mubomegas/2 << endl;
    cout << "matrix.T_6_1()" << "\t" << + 1/sqrt(2)*hbar2Ix/2*sqrt(omegab/omegas)*sqrt(mub/mu) << endl;
    cout << "matrix.T_6_2()" << "\t" << + 1/sqrt(2)*hbar2Ix/2*sqrt(omegas/omegab)*sqrt(mu/mub) << endl;
    cout << "matrix.T_7_1()" << "\t" << hbar2Iz << endl;

    cout << "matrix.Tr()" << "\t" << + B << endl;
    cout << "hbar2Ix" << "\t" << hbar2Ix << endl;
    cout << "matrix.Tr_1()" << "\t" << (hbar2Iz - hbar2Ix) << endl;
    cout << "matrix.Vs(kzz, az)" << "\t" <<kzz << endl;
    cout << "matrix.Vb()" << "\t" <<bhbaromega/4 << endl;
    cout << "V3" << "\t" << V3 << endl;
    cout << "V4" << "\t" << V4 << endl;
    cout << "matrix.Vsb(kxxz, az)" << "\t" << kxxz*Vbcoeff << endl;
    cout << "VR3" << "\t" << VR3*sqrt( Vscoeff ) << endl;
    cout << "VR4" << "\t" << VR4*sqrt( Vscoeff ) << endl;
    cout << "VRR3" << "\t" << VRR3*( Vscoeff ) << endl;
    cout << "VRR4" << "\t" << VRR4*( Vscoeff ) << endl;
    cout << "Vrho3" << "\t" << Vrho3*Vbcoeff << endl;
    cout << "Vrho4" << "\t" << Vrho4*Vbcoeff << endl;
    cout << "Vam31" << "\t" << Vam31*Vbcoeff << endl;
    cout << "VamR31" << "\t" << VamR31*Vbcoeff*pow(Re, -6) << endl;
    cout << "Vam31__" << "\t" << Vam31*sqrt(2*Vbcoeff)/Re << endl;
    cout << "VamR31__" << "\t" << VamR31*sqrt(2*Vbcoeff)/Re*pow(Re, -6) << endl;
    cout << "omegas" << "\t" << omegas << endl;
    cout << "omegab" << "\t" << omegab << endl;
    cout << "Vbcoeff" << "\t" << Vbcoeff << endl;
    cout << "Vscoeff" << "\t" << Vscoeff << endl;

    MatrixParameter dim(m, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
    //vector<vector<double>> result(dim.matrix_size(), vector<double>(dim.matrix_size()));
    Eigen::MatrixXd result = Eigen::MatrixXd::Ones(dim.matrix_size(), dim.matrix_size());
    vector<int> drow(1, 0);
    vector<int> ddrow(1, 0);
    int dddrow= 0;
    for (int mBrow = - j; mBrow <= j; mBrow++) {                                                    // assign row mB
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
            for (int mBcolumn = - j; mBcolumn <= j; mBcolumn++) {                                       // assign column mB
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
                    if (d != dcolumn[dcolumn.size() - 2]) { //dcolumn.size() - 2 > 0の条件を追加したらエラーが消えるかもしれない．
                        ddcolumn.push_back(d);
                        dddcolumn += ddcolumn[ddcolumn.size() - 2];
                    }

                    #pragma omp parallel for num_threads(12)
                    for (int jrow = abs(mBrow); jrow <= j; jrow++) {                                            // assign row j
                        for (int krow = - jrow; krow <= jrow; krow++) {                                         // assign row k
                            for (int n0row = 0; n0row <= n0; n0row++) {                                         // assign row n0
                                for (int jcolumn = abs(mBcolumn); jcolumn <= j; jcolumn++) {                                                // assign column j      
                                    for (int kcolumn = - jcolumn; kcolumn <= jcolumn; kcolumn++) {                                          // assign column k
                                        for (int n0column = 0; n0column <= n0; n0column++) {                                                // assign column n0
                                            // matrix element
                                            MatrixElement matrix(m, jrow, mBrow, krow, n0row, vrow, jcolumn, mBcolumn, kcolumn, n0column, vcolumn);
                                            //cout << "5_1" << ",\t" << - hbar2Ix/2*matrix.T_5_1() << '\n' << endl;
                                            //if (matrix.T_5_2() != 0) {cout << "5_2" << ",\t" << - shbaromega*ReIxsqrt2muomegas/2000*matrix.T_5_2() << ",\t" << n0row << ",\t" << vrow << ",\t" << -mBrow << ",\t" << n0column << ",\t" << vcolumn << ",\t" << -mBcolumn << '\n' << endl;}
                                            result(dddrow + matrix.row(), dddcolumn + matrix.column())
                                             = 
                                             (
                                              - bhbaromega/4*(
                                                + (1 + Re2Ix*mu)*matrix.Tb()    //Tb
                                                + ReIxsqrt2mubomegas*sqrt(mub/mu)*matrix.Tb_1()
                                                + hbar2Ix/2*1/shbaromega*mub/mu*matrix.Tb_2()
                                                )
                                              

                                              - shbaromega/4*matrix.Ts()    //Ts
                                                - hbar2Ix/4*omegas/omegab*mu/mub*matrix.Ts_1()
                                            
                                              
                                              + hbar2Iz*matrix.T_3_1()
                                              - hbar2Ix/2*matrix.T_4_1()
                                              //+ Resqrt2mubhbaromegas/2*matrix.T_4_2()
                                              + hbar2Ix/2*matrix.T_5_1()
                                              + shbaromega*ReIxsqrt2mubomegas/2*matrix.T_5_2()
                                              + 1/sqrt(2)*hbar2Ix/2*sqrt(omegab/omegas)*sqrt(mub/mu)*matrix.T_6_1()
                                              + 1/sqrt(2)*hbar2Ix/2*sqrt(omegas/omegab)*sqrt(mu/mub)*matrix.T_6_2()
                                              + hbar2Iz*matrix.T_7_1()

                                              
                                              + B*matrix.Tr()   //Tr
                                              + hbar2Ix*matrix.Tr() + (hbar2Iz - hbar2Ix)*matrix.Tr_1()
                                              
                                              
                                              + kzz*matrix.Vs(kzz, az)  //Vs
                                              + bhbaromega/4*matrix.Vb()    //Vb  //Vxx termは入れてない。
                                            
                                              //Vr
                                              + V3*( matrix.Vr(3, 2) + matrix.Vr(3, -2) ) 
                                              + V4*( sqrt(14)*matrix.Vr(4, 0) - sqrt(5)*(matrix.Vr(4, 4) + matrix.Vr(4, -4)) )
                                            
                                              //Vsb
                                              + kxxz*bhbaromega/(4*kxx)*matrix.Vsb(kxxz, az)

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
                                              + Vam31*sqrt(2*Vbcoeff)/Re*( (matrix.Vam(3, 2, 1) - matrix.Vam(3, 2, -1)) + (matrix.Vam(3, -2, 1) - matrix.Vam(3, -2, -1)) )
                                              //VamR
                                              + VamR31*sqrt(2*Vbcoeff)/Re*( (matrix.VamR(3, 2, 1) - matrix.VamR(3, 2, -1)) + (matrix.VamR(3, -2, 1) - matrix.VamR(3, -2, -1)) )

                                             );
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
    //PRINT_MAT(result);
    return result;
}