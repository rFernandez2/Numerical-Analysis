#include <iostream>
#include <math.h>
#include <fstream>
#include <string>
#include "gnuplot_i.hpp"

double const MAXLEVEL = 30;
double const MINLEVEL = 5;
double const EXACTVAL = 0.00000000001;

//f(x) = x^2
double f1(double x) {
    return x * x;
}

//f(x) = sqrt(sin(pi * x / 2))
double f2(double x) {
    return sqrt(sin(M_PI * x / 2));
}

//f(x) = x^4
double f3(double x) {
    return pow(x, 4);
}

//f(x) = {1 if x < 0.5, 0 if x >= 0.5}
double f4(double x) {
    return (x < 0.5) ? 1 : 0;
}


double adapTrap(double (*f)(double), double low, double hi, double err, int currLevel,
                double f_lo, double f_hi, int &funcEvals) {
    double currErr;
    double h = (hi+low) / 2;
    double f_h = f(h);
    double coarseApprox = ((f_hi + f_lo) / 2) * (hi-low);
    double fineApprox = (hi - low) / 4 * (f_lo + 2 * f_h + f_hi);
    currErr = fabs(fineApprox - coarseApprox);
    funcEvals ++;
    if((currErr < 3 * err && currLevel >= MINLEVEL) || currLevel >= MAXLEVEL) {
        //printf("%.7f\t%.7f\n", low, -log(fabs(hi-low))); used for plot c by redirecting stdout
        return (4 * fineApprox - coarseApprox) / 3;
    }
    currLevel ++;
    return adapTrap(f, low, h, err / 2, currLevel,f_lo, f_h, funcEvals) +
           adapTrap(f, h, hi, err / 2, currLevel,f_h, f_hi, funcEvals);
}

int main() {
    double tolerances[9] = {.001, .0003, .0001, .00003, .00001, .000003,
                          .000001,.0000003, .0000001};
    double log_evals[4][9];
    double actual_err[4][9];
    double log_tols[9];
    for(int i = 0; i < 9; i ++) {
        log_tols[i] = log(1/tolerances[i]);
    }
    double approxs[4] = {0,0,0,0};
    double exact[4] = {0,0,0,0};
    double lo[4] = {0,0,-1,0};
    double hi[4] = {2,1,1,1};
    double (*f[4]) (double x);
    f[0] = &f1;
    f[1] = &f2;
    f[2] = &f3;
    f[3] = &f4;

    for(int i = 0; i < 4; i ++) {
        int funcEvals = 0;
        exact[i] = adapTrap(f[i], lo[i], hi[i], EXACTVAL, 1, f[i](lo[i]), f[i](hi[i]), funcEvals);
        for(int j = 0; j < 9; j ++) {
            funcEvals = 0;
            approxs[i] = adapTrap(f[i], lo[i], hi[i], tolerances[j], 1, f[i](lo[i]), f[i](hi[i]), funcEvals);
            log_evals[i][j] = log(funcEvals);
            actual_err[i][j] = fabs(exact[i] - approxs[i]);
            printf("Function #%d, Func Evals: %d, Approximation: %.5f, \"Exact\": %.5f, Difference: %.5f\n",
                   i + 1, funcEvals, approxs[i], exact[i], actual_err[j][i]);
        }
    }

    //Below code used to make plots. Gnuplot is hard so some was done manually (i.e. plots c)
    /*
    Gnuplot g1("lines");
    g1.set_title("Function 4 Plot A");
    std::vector<double> x, y;
    for (int j = 0; j < 9; j++) {
        x.push_back(log_tols[j]);
        y.push_back(log_evals[3][j]);
    }
    g1.set_smooth().plot_xy(x, y, "logvals vs. logtols"); */

    /*Gnuplot g2("lines");
    g2.set_title("Function 1 Plot B (4*fine - coarse) / 3");
    std::vector<double> x,y;
    for(int i = 0; i < 9; i ++) {
        x.push_back(log_tols[i]);
        y.push_back(log(1 / actual_err[0][i]));
        //y.push_back(0); used for function 4
    }
    g2.set_smooth().plot_xy(x, y, "logvals vs. logerror");*/

/*
    Gnuplot g3("lines");
    g3.set_title("Function 1 Plot C");
    std::vector<double> x,y;
    for(int i = 0; i < 9; i ++) {
        x.push_back(log_tols[i]);
        //y.push_back(log(1 / actual_err[3][i]));
        y.push_back(0);
    }
    g2.set_smooth().set_yrange(-0.01, 0.01).plot_xy(x, y, "logvals vs. logerror");*/

    return 0;
}
