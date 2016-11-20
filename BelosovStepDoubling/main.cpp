 #include <iostream>
#include <math.h>

#include <fstream>
#include "stiff.h"

 std::ofstream tracef1("trace1.0e-2");


 Matrix f( double t, const Matrix& y, Matrix& fa, Matrix& df ) {
     Matrix ans(y.row, y.col);
     ans(0,0) = -y(0,0);

     fa(0,0) = y(0,0);

     df(0,0) = -1;
     return ans;
 }

 double fnorm(const Matrix& y) {
     return fabs(y(0,0));
 }

Matrix fBel( double t, const Matrix& y, Matrix& fa , Matrix& df  ) {
    Matrix ans(y.row, y.col);
    ans(0,0) = 77.27 * (y(1,0) - y(1,0) * y(0,0) + y(0,0) - 8.375e-6 * pow(y(0,0), 2));
    ans(1,0) = (1 / 77.27) * (-y(1,0) - y(1,0)*y(0,0) + y(2,0));
    ans(2,0) = 0.161 * (y(0,0) - y(2,0));

    df(0,0) = 77.27 * (-1 * y(1,0) + 1 - 2 * 8.375e-6 * y(0,0)); df(0,1) = 77.27 * (1 - y(0,0)); df(0,2) = 0;
    df(1,0) = -1 / 77.27 * y(1,0); df(1,1) = -1 / 77.27 * (1 + y(0,0)); df(1,2) = 1 / 77.27;
    df(2,0) = 0.161; df(2,1) = 0; df(2,2) = -0.161;

    fa(0,0) = 77.27 * (fabs(y(1,0)) + fabs(y(1,0)) * fabs(y(0,0)) + fabs(y(0,0)) + (8.375 * pow(10, -6)) * pow(y(0,0), 2));
    fa(1,0) = (1 / 77.27) * (fabs(y(1,0)) + fabs(y(1,0)) * fabs(y(0,0)) + fabs(y(2,0)));
    fa(2,0) = 0.161 * (fabs(y(0,0)) + fabs(y(2,0)));

    return ans;
}

double belnorm(const Matrix& y) {
    return fabs(y(0,0)) / (1.25 * pow(10,5)) + fabs(y(1,0)) / 1800 + fabs(y(2,0)) / (3 * pow(10,4));
}

void testBelosov() {
    Matrix Ynew(3,1);
    Ynew(0,0) = 4; Ynew(1,0) = 1.1;   Ynew(2,0) = 4;
    Matrix Yold(Ynew);
    simTime stime; stime.endTime = 700;
    double tols[3] = {1.0e-3, 1.0e-4, 1.0e-5};
    for(int i = 0; i < 6; i ++) {
        //if(i % 2 == 0)
            //RESIDUAL_LIM = 200;
        stime.tol = tols[i];
        advance(stime, Ynew, Yold, fBel, belnorm);
    }

}

 void testExp() {
     Matrix Ynew(1,1);
     Ynew(0,0) = 2;
     Matrix Yold(Ynew);
     simTime stime; stime.endTime = 2;
     advance(stime, Ynew, Yold, f, fnorm);
 }


int main() {
    //testExp();
    testBelosov();
    return 0;
}
