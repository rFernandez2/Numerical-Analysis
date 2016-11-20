#include "qrsolve.h"

double ERROR = 1e-4;

double eNorm(const Matrix &A) {
    double temp = 0;
    for(uint i = 0; i < A.row; i ++) {
        temp += pow(A(i,0),2);
    }
    return sqrt(temp);
}

void factor(Matrix &A, Matrix &Q, Matrix &R) {
    /*Matrix e0(A.row + 1, 1); e0(0,0) = 1; //first standard basis vector R^(n+1)
    Matrix V(transpose(A) - eNorm(A) * e0);
    */
    /*double mag, beta;
    uint m = A.row; uint n = A.col;
    Matrix u(m, 1), v(m, 1);
    Matrix P(m, m), I(m, m);
    Q = Matrix(m, m);
    for(uint i = 0; i < I.row; i ++) {
        I(i, i) = 1;
        P(i, i) = 1;
        Q(i, i) = 1;
    }

    R = A;

    for (uint i = 0; i < n; i++) {
        zero(u);
        zero(v);
        mag = 0.0;
        for (uint j = i; j < m; j++) {
            u(j, 0) = R(j, i);
            mag += u(j,0) * u(j,0);
        }
        mag = sqrt(mag);

        beta = u(i,0) < 0 ? mag : -mag;
        mag = 0.0;
        for (uint j = i; j < m; j++) {
            v(j,0) = j == i ? u(j,0) + beta : u(j,0);
            mag += v(j,0) * v(j,0);
        }
        mag = sqrt(mag);

        if(mag < ERROR) continue;

        for (uint j = i; j < m; j++) v(j,0) /= mag;
        Matrix vt(transpose(v));
        P = I - 2.0 * (v * vt);
        R = P * R;
        Q = Q * P;
    }*/
    /*double s;
    uint m = A.row; uint n = A.col;
    for(uint i = 0; i < n; i ++) {
       for(uint j = 0; j < i; j ++) {
           s = 0;
           for(uint k = 0; k < m; k ++) {
               s += A(k,j) * A(k, i);
               R(j, i) = s;
           }
       }
        for(uint j = 0; j < i - 1; j ++) {
            for(uint k = 0; k < m; k ++) {
                A(k, i) -= A(k, j) * R(j,i);
            }
        }
        s = 0;
        for(uint j = 0; j < m; j ++) {
            s += pow(A(j, i), 2);
        }
        R(i,i) = sqrt(s);
        for(uint j = 0; j < m; j ++) {
            A(j, i) /= R(i, i);
        }
    }*/
    double s;
    uint m = A.row; uint n = A.col;
    for(uint i = 0; i < n; i ++) {
        s = 0;
        for(uint j = 0; j < m; j ++) {
            s += pow(A(j, i), 2);
        } R(i,i) = sqrt(s);
        for(uint j = 0; j < m; j ++) Q(j, i) = A(j, i) / R(i,i);
        for(uint j = i + 1; j < n; j ++) {
            s = 0;
            for(uint k = 0; k < m; k ++) s+= A(k, j) * Q(k, i);
            R(i, j) = s;
            for(uint k = 0; k < m; k ++) A(k, j) -= R(i, j) * Q(k, i);
        }
    }
}