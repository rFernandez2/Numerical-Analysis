// interface for tools for a Galerkin method for 2 point BVPs
// that uses continuous piecewise linears

#include "matrix.h"
#include <iostream>
#include <vector>

typedef double (*d2dfp)(double);

enum BCtype {Dirichlet, Neumann};

// compute the Galerkin approximation of the 2 point BVP
Matrix galerkinBVP( d2dfp k,   // diff.eq. is 
                    d2dfp v,   // -(k y')' + v y' + g y = f
                    d2dfp g,
                    d2dfp f,
                    const std::vector<double>& x, // partition of [a.b]
                    double ya,  // boundary data and types for a & b
                    BCtype bca,
                    double yb,
                    BCtype bcb );

// compute the tridiag and the RHS for NN problem
tridiag galerkinMats(d2dfp k,   // diff.eq. is 
                     d2dfp v,   // -(k y')' + v y' + g y = f
                     d2dfp g,
                     d2dfp f,
                     const std::vector<double>& x, // partition of [a.b]
                     Matrix& rhs);

// write out the x,y pairs that define the solution
void writexy(const std::vector<double>& x,
             const Matrix& y, std::ostream & out);

// check vector -- it should be monotone strictly monotone
bool strictlyMonotone( std::vector<double> x );
