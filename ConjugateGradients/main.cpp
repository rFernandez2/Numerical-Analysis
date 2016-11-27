// test CG

#include "conjugateGrad.h"
#include <cmath>
#include <iostream>
#include <fstream>

double sinsin(double x, double y)
{
  return sin(M_PI*2*x)*sin(M_PI*3*y);
}

double bumpbump(double x, double y)
{
  return -2*(x*(1-x)+y*(1-y));
}

double zero(double x, double y)
{
  return 0.0;
}

int main() {
    int N = 10;
    Matrix u(N, N); // u is zero -- initial guess at the solution
    Matrix r(N, N);
    double dx = 1.0 / (N - 1);
    fill(r, sinsin, dx); // fill r with sin(pi*2*x)*sin(pi*3*y)
    Matrix s(r);
    double normRessq = (r, r);

    std::cout << "Test case 1: The initial norm squared of the residual is "
              << normRessq << std::endl;
    stepCG(u, s, r, normRessq);
    std::cout << "After one step of CG the norm squared of the residual is "
              << normRessq << std::endl;
    stepCG(u, s, r, normRessq);

    u = 0.0 * r; // initial guess is zero
    fill(r, bumpbump, dx);
    s = r;
    normRessq = (r, r);
    std::cout << "Test case 2: The initial norm squared of the residualis "
              << normRessq << std::endl;
    int max_iter = 10;
    for (int itr = 1; itr <= max_iter; itr++) {
        stepCG(u, s, r, normRessq);
        std::cout << " after interation " << itr
                  << " the normsquared of the residual is " << normRessq
                  << std::endl;
    }
    std::ofstream out("ex2");
    for(uint i = 0; i < u.Row(); i ++) {
        for(uint j = 0; j < u.Col(); j ++) {
            out << i / 9.0 << '\t' << j / 9.0 << '\t' << u(i,j) << '\n';
        }
        out << '\n';
    }

    out.close();

    system("gnuplot ./plot.gnu");

    return 0;
}
