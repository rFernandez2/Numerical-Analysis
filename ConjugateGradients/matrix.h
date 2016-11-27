#ifndef MATRIXH
#define MATRIXH
#include <cmath>
#include <iostream>

#include "utility.h"

class Matrix{
  // This is a row x col full matrix
 int row, col;
 double *m;

public:
 Matrix(int r , int c=1);       // constructors & destructor
 Matrix( const Matrix& );
 ~Matrix(){delete [] m;};

 int Row() const  {return row;}
 int Col() const {return col;}
 inline double& operator()(int i, int j=0)
 { if( 0<=i && i<row && 0<=j && j<col)
     return (m[i*col+j]);
   else error("Matrix[]: index bust");
   return m[0]; //get rid of a compiler warning
 }

 inline double operator()(int i, int j=0)const
 { if( 0<=i && i<row && 0<=j && j<col)
     return (m[i*col+j]);
   else error("Matrix[]: index bust");
   return m[0]; //get rid of a compiler warning
 }

 Matrix& operator=( const Matrix& B);   // *this=B
 Matrix& operator+=( const Matrix& B);   // *this+=B
 Matrix& operator-=( const Matrix& B);   // *this-=B


 friend Matrix operator+(const Matrix& A, const Matrix& B); // A+B
 friend Matrix operator*(const Matrix& A, const Matrix& B); // A*B
 friend Matrix operator*( double a, const Matrix& B); // a*B
 friend Matrix operator-(const Matrix& A, const Matrix& B); // A-B
 friend Matrix transpose(const Matrix&); // return A^T
 friend Matrix abs(const Matrix&); // take abs of entries
 friend double operator,(const Matrix& A, const Matrix& B); // A:B
};

std::ostream& operator<<(std::ostream& os, const Matrix& A);


#endif
