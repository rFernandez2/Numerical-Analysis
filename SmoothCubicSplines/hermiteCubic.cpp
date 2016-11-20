//
// Created by Roberto Fernandez on 10/23/16.
//

#include "hermiteCubic.h"

double hermiteCubic::evaluate(double x) {
    double temp; //stores f, f', f''
    int l = (int)(lower_bound(p, p + n, x) - p); //finds left endpoint of partition where x is
    int b = l - 1; //left endpoint of interval right before x
    double dp0 = p[l] - p[b];
    double x0 = (x-p[b]) / dp0;
    temp = v[b] * v0(x0)  + d[b] * s0(x0) * dp0
              + v[l] * v1(x0) + d[l] * s1(x0) * dp0; //Calculates g(x)
    return temp;
}

std::vector<double> smoothCubic_DEC(hermiteCubic h) {
    std::vector<double> temp;
    TridiagonalMatrix t(h.n);
    t(0,0) = 1; t(0,1) = 0; t(h.n-1, h.n - 2) = 0; t(h.n - 1, h.n - 1) = 1;

    double dt0, dt1;

    for(uint i = 1; i < h.n - 1; i ++) {//calculates the A matrix to be solved
        dt0 = h.p[i] - h.p[i - 1];
        dt1 = h.p[i + 1] - h.p[i];
        t(i, i - 1) = dt1;
        t(i,i) = 2 * (dt0 + dt1);
        t(i, i + 1) = dt0;
    }

    Matrix d(h.n, 1);
    Matrix y(h.n, 1);

    dt0 = h.p[1] - h.p[0];
    double dv0, dv1;
    dv0 = h.v[1]-h.v[0];


    y(0,0) = h.d[0];
    for(uint i = 1; i < h.n - 1; i ++) { //sets up y matrix to solve Ax=y
        dt1 = h.p[i+1] - h.p[i];
        dv1 = h.v[i + 1] - h.v[i];
        y(i,0) = 3.0 / (dt0 * dt1) * (pow(dt0, 2) * dv1 + pow(dt1,2) * dv0);
        dt0 = dt1;
        dv0 = dv1;
    }
    y(h.n-1, 0) = h.d[h.n-1];

    //Solves
    factor(t);
    d = solve(t, y);

    //Copies over result from resultant matrix
    for(uint i = 0; i < h.n; i ++) {
       temp.push_back(d(i,0));
        h.d[i] = d(i,0);
    }
    return temp;
}

std::vector<double> not_a_knot(hermiteCubic h) {
    if(h.n < 4)
        error("Must have at least 4 points to do not_a_knot.");

    uint n = h.n - 1; //Easier indexing
    std::vector<double> temp;
    TridiagonalMatrix t(h.n, h.n);
    double dt0, dt1;

    dt0 = h.p[1] - h.p[0];
    dt1 = h.p[2] - h.p[1];
    t(0,0) = dt1; //Writes C^3 at 1 condition
    t(0,1) = dt1 + dt0;

    for(uint i = 1; i < n; i ++) { //sets up A matrix
        dt0 = h.p[i] - h.p[i - 1];
        dt1 = h.p[i + 1] - h.p[i];
        t(i, i - 1) = dt1;
        t(i,i) = 2 * (dt0 + dt1);
        t(i, i + 1) = dt0;
    }
    t(n, n-1) = h.p[n] - h.p[n-2]; //Writes C^3 at n-1 condition
    t(n, n) = h.p[n-1] - h.p[n-2];

    Matrix d(h.n, 1);
    Matrix y(h.n, 1);
    dt0 = h.p[1] - h.p[0];
    dt1 = h.p[2] - h.p[1];
    double dv0, dv1;
    dv0 = h.v[1]-h.v[0];
    double num = (dt0 + 2 * (dt0 + dt1)) * dt1 * (h.v[1] - h.v[0]) / dt0 + (pow(dt0,2)) * (h.v[2] - h.v[1]) / dt1;
    y(0,0) = num / (dt0 + dt1); //Calculates C^3 at 1 condition (sorry about mess!)

    for(uint i = 1; i < h.n - 1; i ++) { //Interior of y matrix
        dt1 = h.p[i+1] - h.p[i];
        dv1 = h.v[i + 1] - h.v[i];
        y(i,0) = 3.0 / (dt0 * dt1) * (pow(dt0, 2) * dv1 + pow(dt1,2) * dv0);
        dt0 = dt1;
        dv0 = dv1;
    }

    dt1 = h.p[n] - h.p[n-1];
    dt0 = h.p[n-1] - h.p[n-2];
    num = (pow(dt1,2)) * (h.v[n-1] - h.v[n-2]) / dt0 + (2 * (dt0 + dt1) + dt1) * dt0 * (h.v[n] - h.v[n-1]) / dt1;
    y(n,0) = num / (dt0 + dt1); //C^3 at n-1 condition


    factor(t);
    d = solve(t, y);

    for(uint i = 0; i < h.n; i ++) {
        temp.push_back(d(i,0));
        h.d[i] = d(i,0);
    }
    return temp;
}

