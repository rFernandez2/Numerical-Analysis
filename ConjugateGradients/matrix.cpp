#include "matrix.h"

Matrix::Matrix(int r, int c){ // set up a matrix of zeros
   if (r<=0 || c<=0)error("bad matrix indeces \n");
   row=r;
   col=c;
   m=new double[r*c];
   for(int i=0; i<r*c; i++) m[i]=0.0;
 };

Matrix::Matrix( const Matrix& A){ // build a matrix from another one
   row = A.row;
   col = A.col;
   m = new double[row*col];
   for(int i=0; i< row*col; i++) m[i] =A.m[i];
 };

Matrix& Matrix::operator=( const Matrix& B){   // *this=B
   // note that this does not follow the standard form for
   // writing assignment operators. Instead it says it is
   // an error to write A = B unless A & B have the same dimensions.

   if( row!=B.row || col!=B.col) error("Matrix=: incompatible sizes");
   int i, rc = row*col;
   for( i=0; i<rc; i++ ) m[i] = B.m[i];
   return *this;
 }

Matrix& Matrix::operator+=( const Matrix& B){   // *this+=B
   if( row!=B.row || col!=B.col) error("Matrix+=: incompatible sizes");
   int i, rc = row*col;
   for( i=0; i<rc; i++ ) m[i] += B.m[i];
   return *this;
 }

Matrix& Matrix::operator-=( const Matrix& B){   // *this-=B
   if( row!=B.row || col!=B.col) error("Matrix-=: incompatible sizes");
   int i, rc = row*col;
   for( i=0; i<rc; i++ ) m[i] -= B.m[i];
   return *this;
 }

Matrix operator+(const Matrix& A, const Matrix& B){ // add two Matrices
 if (A.row!=B.row || A.col!=B.col)error("Matrix+: incomatable sizes \n");
 Matrix C(A.row, A.col);
 int i,rc=A.row*A.col;
 for(i=0; i<rc; i++)
    C.m[i] = A.m[i] + B.m[i];
 return C;
}

Matrix operator*(const Matrix& A, const Matrix& B){
 // Matrix multiplication
 if (A.col!=B.row)error("Matrix*: incomatable sizes \n");
 Matrix C(A.row, B.col);
 int i,j,k;
 for(i=0; i<A.row; i++)
    for(j=0; j<B.col; j++)
       for(k=0; k<A.col; k++)
          C.m[i*C.col+j]+=A.m[i*A.col + k] * B.m[k*B.col+j];
 return C;
}

Matrix operator*( double a, const Matrix& B){
 // a*B
 Matrix C( B);
 int i,s=C.row*C.col;
 for(i=0; i<s; i++)
   C.m[i] *=a;
 return C;
}

Matrix operator-(const Matrix& A, const Matrix& B){
 // Matrix subtraction
 if (A.row!=B.row || A.col!=B.col)error("Matrix-: incomatable sizes \n");
 Matrix C(A.row, A.col);
 int i,rc=A.row*A.col;
 for(i=0; i<rc; i++)
    C.m[i]=A.m[i] - B.m[i];
 return C;
}

std::ostream& operator<<(std::ostream& os, const Matrix& A){
 // output for small matrices
 int i, j, row=A.Row(), col=A.Col();
 for(i=0; i<row; i++){
    for(j=0; j<col; j++)
       os << A(i,j) << "\t";
    os << "\n";
   }
 return os;
}

Matrix transpose (const Matrix& A){
 Matrix C(A.col, A.row);
 int i,j;
 for(i=0; i<A.row; i++)
    for(j=0; j<A.col; j++)
       C.m[j*C.col + i]=A.m[i*A.col + j];
 return C;
}

Matrix abs (const Matrix& A){
  Matrix C(A);
 int i,k=A.col*A.row;
 for(i=0; i<k; i++)
   C.m[i] = abs(C.m[i]);
 return C;
}

// A:B
double operator,(const Matrix& A, const Matrix& B){
  // sanity check
  if( A.row!=B.row || A.col!=B.col )
    error("operator,: shape mismatch");

  double s=0.0;
  int j, n = A.row*A.col;
  for( j=0; j<n; j++)
    s += A.m[j]*B.m[j];

  return s;
}
