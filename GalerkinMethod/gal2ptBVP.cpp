//
// Created by Roberto Fernandez on 10/30/16.
//
#include "gal2ptBVP.h"

tridiag galerkinMats(d2dfp k,   // diff.eq. is
                     d2dfp v,   // -(k y')' + v y' + g y = f
                     d2dfp g,
                     d2dfp f,
                     const std::vector<double>& x, // partition of [a.b]
                     Matrix& b,
                     double ya,
                     double yb) {
    uint l = (uint) x.size();
    TridiagonalMatrix A(l, 0, 0, 0);
    double dx, xm; // dx = x_i+1-x_i and xm = midpoint of current interval [x_i-1,x_i]
    double km, vm, gm, fm; //vals at xm for each given function
    for(uint i = 1; i < l; i ++) {
        //solving sub-problem of 2x2 matrix and inserting into main tridiagonal and column matrices
        dx = x[i] - x[i - 1];
        xm = x[i - 1] + dx / 2;
        km = k(xm); vm = v(xm); gm = g(xm); fm = f(xm);
        A(i - 1,i - 1) += km / dx - vm / 2 + gm * dx / 3;
        A(i - 1, i) += -km / dx + vm / 2 + gm * dx / 6;
        A(i, i - 1) += -km / dx - vm / 2 + gm * dx / 6;
        A(i, i) += km / dx + vm / 2 + gm * dx / 3;

        b(i - 1, 0) += fm * dx / 2;
        b(i, 0) += fm * dx / 2;
    }
    b(0,0) -= k(x[0]) * ya; //Includes left Neumann condition
    b(l-1, 0) += k(x[l - 1]) * yb; //Includes right Neumann condition
    return A;
}

// compute the Galerkin approximation of the 2 point BVP
Matrix galerkinBVP( d2dfp k,   // diff.eq. is
                    d2dfp v,   // -(k y')' + v y' + g y = f
                    d2dfp g,
                    d2dfp f,
                    const std::vector<double>& x, // partition of [a.b]
                    double ya,  // boundary data and types for a & b
                    BCtype bca,
                    double yb,
                    BCtype bcb ) {
    uint l = (uint) x.size();
    Matrix b(l,1);
    Matrix y(l,1);
    TridiagonalMatrix A(galerkinMats(k, v, g, f, x, b, ya, yb)); //Setting up Ay=b

    if(bca == Dirichlet) { //DN or DD conditions
        A(0,0) = 1; A(0, 1) = 0; //Setting y(a)=y_a
        b(0,0) = ya;
    }
    if(bcb == Dirichlet) { //ND or DD conditions
        A(l - 1, l - 2) = 0; A(l - 1, l - 1) = 1; //Setting y(b)=y_b
        b(l - 1, 0) = yb;
    }

    factor(A);
    y = solve(A, b);

    return y;
}

void writexy(const std::vector<double>& x, const Matrix& y, std::ostream & out) {
    uint l = (uint) x.size();
    for(uint i = 0; i < l; i ++) {
        out << x[i] << '\t' << y(i,0) << std::endl;
    }
}
