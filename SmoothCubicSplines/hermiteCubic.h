//
// Created by Roberto Fernandez on 10/23/16.
//

#include "matrix.h"
#include <vector>
#include <math.h>

#ifndef SMOOTHCUBICSPLINES_HERMITECUBIC_H
#define SMOOTHCUBICSPLINES_HERMITECUBIC_H

using namespace std;

static inline double v0(double x) {
    return 1 - 3 * pow(x,2) + 2 * pow(x,3);
}

static inline double v1(double x) {
    return 3 * pow(x,2) - 2 * pow(x,3);
}

static inline double s0(double x) {
    return x * pow(x - 1, 2);
}

static inline double s1(double x) {
    return pow(x,2) * (x - 1);
}

class hermiteCubic {
public:
    uint n;
    double *d, *p, *v; //d = g' at points p. p = interval endpoints, v = g at points pi

    double evaluate(double x);
    hermiteCubic(uint num, double *der, double *par, double *val) {
        n = num;
        p = par;
        v = val;
        d = der;
    }
};

vector<double> smoothCubic_DEC(hermiteCubic f);

std::vector<double> not_a_knot(hermiteCubic h);


#endif //SMOOTHCUBICSPLINES_HERMITECUBIC_H
