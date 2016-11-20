#ifndef MATRIXH
#define MATRIXH
#include <iostream>

#include "utility.h"

typedef unsigned int uint;

class Matrix{
  // This is a row x col matrix
    public: uint row, col;
 double *m;

public:
 Matrix(uint r , uint c=1);       // constructors & destructor
 Matrix( const Matrix& );
 ~Matrix(){delete [] m;};

 uint Row() const  {return row;}
 uint Col() const {return col;}
 inline double& operator()(uint i, uint j=0)
 { return (m[i*col+j]);}
 inline double operator()(uint i, uint j=0)const
 { return (m[i*col+j]);}

 // error if size mismatch
 Matrix& operator=( const Matrix& B);   // *this=B

 friend Matrix operator+(const Matrix& A, const Matrix& B); // A+B
 friend Matrix operator*(const Matrix& A, const Matrix& B); // A*B
 friend Matrix operator*( double a, const Matrix& B);       // a*B
 friend Matrix operator-(const Matrix& A, const Matrix& B); // A-B
 friend void factor(Matrix& A); // Factor the square matric A in place
  // solve Ax=y using its factored form
 friend Matrix solve(const Matrix& AF, const Matrix& y);
};

std::ostream& operator<<(std::ostream& os, const Matrix& A);

// square tridiagonal matrices
class TridiagonalMatrix;
typedef TridiagonalMatrix tridiag;
class TridiagonalMatrix{
  uint row;
  double *b,*d,*a; // mnemonic below, diag, above
 public:
  // build tridiag with constant diagonals
  TridiagonalMatrix( uint r, double bv=0.0, double dv=0.0, double av=0.0 );
  TridiagonalMatrix( const tridiag & A );
  ~TridiagonalMatrix(){ delete [] b;};

  uint Row() const  {return row;}
  uint Col() const {return row;} // it is square
  double& operator()(uint i, uint j);
  double operator()(uint i, uint j) const;

 // error if size mismatch
  TridiagonalMatrix& operator=( const TridiagonalMatrix& B);

  friend tridiag operator+(const tridiag& A, const tridiag& B);
  friend Matrix operator*(const tridiag& A, const Matrix& B);
  friend tridiag operator*( double a, const tridiag& B);
  friend tridiag operator-(const tridiag& A, const tridiag& B);
  friend void factor( tridiag& A ); // factor in place
  // solve using the factors
  friend Matrix solve(const tridiag& AF, const Matrix& y);

};

#endif






