// interface for CG on a 5-point laplace operator

// Here the vectors involved in the linear algebra are
// doubly indexed because they represent 2-d functions.
// So U(i,j) is the value of our approximation to u(i*dx,j*dy).

#include "matrix.h"

// multiply the vector u by the five point laplacian
// The problem has homogenous Dirichle boundary values
// so the matrix C has rows that look like the identity
// on the boundary. The vector u should be zero on the
// boundary.
// It is assumed that the mesh is on (0,1)x(0,1) and
// that the mesh is NxN
void mul5ptLaplace( const Matrix& u, Matrix& Cu );

// Take one step of CG
void stepCG( Matrix& y, Matrix& s, Matrix& r, double& normResSq );

// fill a vector with function values
typedef double (*pt2dbl)( double x, double y);
void fill( Matrix& u, pt2dbl f, double dx);
