#include "matrix.h"

Matrix::Matrix(uint r, uint c){ // set up a matrix of zeros
   if (r==0 || c==0)error("bad matrix indeces \n");
   row=r;
   col=c;
   m=new double[r*c];
   for(uint i=0; i<r*c; i++) m[i]=0.0;
 };

Matrix::Matrix( const Matrix& A){ // build a matrix from another one
   row = A.row;
   col = A.col;
   m = new double[row*col];
   for(uint i=0; i< row*col; i++) m[i] =A.m[i];
 };

Matrix& Matrix::operator=( const Matrix& B){   // *this=B
   // note that this does not follow the standard form for
   // writing assignment operators. Instead it says it is
   // an error to write A = B unless A & B have the same dimensions.

   if( row!=B.row || col!=B.col) error("Matrix=: incompatible sizes");
   uint i, rc = row*col;
   for( i=0; i<rc; i++ ) m[i] = B.m[i];
   return *this;
 }

Matrix operator+(const Matrix& A, const Matrix& B){ // add two Matrices
 if (A.row!=B.row || A.col!=B.col)error("Matrix+: incomatable sizes \n");
 Matrix C(A.row, A.col);
 uint i,rc=A.row*A.col;
 for(i=0; i<rc; i++)
    C.m[i] = A.m[i] + B.m[i];
 return C;
}


Matrix operator*(const Matrix& A, const Matrix& B){
 // Matrix multiplication
 if (A.col!=B.row)error("Matrix*: incomatable sizes \n");
 Matrix C(A.row, B.col);
 uint i,j,k;
 for(i=0; i<A.row; i++)
    for(j=0; j<B.col; j++)
       for(k=0; k<A.col; k++)
         C(i,j) +=A(i,k) * B(k,j);
 return C;
}

Matrix operator*( double a, const Matrix& B){
 // a*B
 Matrix C( B);
 uint i,s=C.row*C.col;
 for(i=0; i<s; i++)
   C.m[i] *=a;
 return C;
}

Matrix operator-(const Matrix& A, const Matrix& B){
 // Matrix subtraction
 if (A.row!=B.row || A.col!=B.col)error("Matrix-: incomatable sizes \n");
 Matrix C(A.row, A.col);
 uint i,rc=A.row*A.col;
 for(i=0; i<rc; i++)
    C.m[i]=A.m[i] - B.m[i];
 return C;
}

std::ostream& operator<<(std::ostream& os, const Matrix& A){
 // output for small matrices
 uint i, j, row=A.Row(), col=A.Col();
 for(i=0; i<row; i++){
    for(j=0; j<col; j++)
       os << A(i,j) << "\t";
    os << "\n";
   }
 return os;
}

void factor( Matrix& A ){

  // To solve A x = B we first factor A in place
  // It is assumed that A can be factored without pivoting.
  // The factors are saved in A.

 if( A.row!=A.col )
    error("Matrix factor: incompatible sizes\n");

    double temp;
    for(int diag = 0; diag < A.row; diag ++) {
        temp = 1.0 / A(diag,diag);
        A(diag,diag) = temp; //reciprocal of diagonal entry
        for(int r = diag + 1; r < A.row; r ++) { //multiplies current row with reciprocal of A(i,i)
            A(diag,r) *= temp;
        }
        for(int j = diag + 1; j < A.row; j ++) { //multiplies rest of rows, j traverses rows
            for(int k = diag + 1; k < A.col; k ++) { //k traverses columns
                A(j,k) -= A(j,diag) * A(diag, k);
            }
        }
    }

}

Matrix solve(const Matrix& A, const Matrix& B){
  // solve A x = B and return the answer.
  // It is assumed that A has been factored without pivoting,
  // using the function factor() that saves the factors of A
  // in place.

 if( A.row!=A.col || A.row!=B.row )
    error("Matrix solve: incompatible sizes\n");

    Matrix solution(A.row, B.col);
    double temp;
 // Add code here
    for(int col = 0; col < B.col; col ++) { //mirrors operations used to factor A and does them on B
        for (int i = 0; i < A.row; i++) {
            solution(i, col) = B(i, col); //copies over val from B matrix
            for (int j = 0; j < i; j++) {
                solution(i, col) -= A(i, j) * solution(j, col); //Mirrors multiplication done when factoring A
            }
            solution(i, col) *= A(i, i); //Multiplies by row's respective diagonal entry
        }

        //Backsolving. Uses calculated values of x_i to simplify and find current x
        for (int i = A.row - 1; i >= 0; i--) {
            for (int j = A.col - 1; j > i; j--) {
                solution(i, col) -= A(i, j) * solution(j, col);
            }
        }
    }
    return solution;
}

TridiagonalMatrix::TridiagonalMatrix
( uint r, double bv, double dv, double av )
{
  row = r;
  b = new double[3*row-2];
  d = b + row-1;
  a = d + row;
  uint i;
  for( i=0; ; i++ ){
    d[i] = dv;
    if( i+1 == r ) break;
    a[i] = av;
    b[i] = bv;
  }
}

TridiagonalMatrix::TridiagonalMatrix( const tridiag & A )
{
  row = A.row;
  b = new double[3*row-2];
  d = b + row-1;
  a = d + row;
  uint j;
  for(j=0; j<3*row-2; j++ )
    b[j] = A.b[j];
}

double& TridiagonalMatrix::operator()(uint i, uint j)
{
  if(  i>=row || j>= row )
    error( "TridiagonalMatrix: index out of range");
  if( j==i-1 )return b[j];
  if( i==j ) return d[j];
  if( j==i+1 )return a[i];
  error( "TridiagonalMatrix: bad index");
  return d[0]; // get rid of compiler warning
}

double TridiagonalMatrix::operator()(uint i, uint j)const
{
  if(  i>=row || j>= row )
    error( "TridiagonalMatrix: index out of range");
  if( j==i-1 )return b[j];
  if( i==j ) return d[j];
  if( j==i+1 )return a[i];
  error( "TridiagonalMatrix: bad index");
  return d[0]; // get rid of compiler warning
}

TridiagonalMatrix& 
TridiagonalMatrix::operator=( const TridiagonalMatrix& B)
{  // *this=B
  if( row != B.row  )
    error("TridiagonalMatrix operator=: size mismatch");
  uint j;
  for( j=0; j<3*row-2; j++ )
    b[j] = B.b[j];
  return *this;
}

tridiag operator+(const tridiag& A, const tridiag& B)
{ // A+B
  if( A.row != B.row )
    error("Tridiagonal matrix operator+: size mismatch");
  tridiag C(A);
  uint j;
  for( j=0; j<3*A.row-2; j++)
    C.b[j] += B.b[j];
  return C;
}

Matrix operator*(const tridiag& A, const Matrix& B)
{ // A*B
  if( A.row != B.Row() )
    error("tridiag*matrix: size mismatch");
  uint j,k;
  Matrix C(B);
  for( k=0; k<B.Col(); k++ ){
    C(0,k) = A.d[0]*B(0,k);
      for(j=0; j+1<A.row; j++){
        C(j,k) += A.a[j]*B(j+1,k);
        C(j+1,k) = A.d[j+1]*B(j+1,k) + A.b[j]*B(j,k);
      }
  }
  return C;
}

tridiag operator*( double a, const tridiag& B)
{ // a*B
  tridiag C(B);
  uint j;
  for( j=0; j<3*B.row-2; j++ )
    C.b[j] *= a;
  return C;
}

tridiag operator-(const tridiag& A, const tridiag& B)
{ // A-B
  if( A.row != B.row )
    error("Tridiagonal matrix operator+: size mismatch");
  tridiag C(A);
  uint j;
  for( j=0; j<3*A.row-2; j++)
    C.b[j] -= B.b[j];
  return C;
}

void factor( tridiag& A )
{ // Factor A in place using Gaussian
  // elimination without pivoting
    double temp;
    int diag = 0;
    for(; diag < A.row - 1; diag ++) {
        temp = 1.0 / A(diag,diag);
        A(diag,diag) = temp; //reciprocal of diagonal entry
        A(diag, diag + 1) *= temp; //multiplies current row with reciprocal of A(i,i)

        A(diag + 1, diag + 1) -= A(diag + 1, diag) * A(diag, diag + 1); //multiples next diagonal
    }
    A(diag, diag) = 1.0 / A(diag, diag); //does last element (1x1 matrix) since diag = A.row - 1
}

Matrix solve(const tridiag& A, const Matrix& y)
{ // solve Ax=y
  // It is assumed that A has been factored by the
  // function factor() and that the factors are in A.

    Matrix solution(A.row, y.col);
    double temp;
    // Add code here
    for(int col = 0; col < y.col; col ++) { //mirrors operations used to factor A and does them on B
        solution(0, col) = y(0, col) * A(0,0); //does the initial case since no operations other than multiply
        for (int i = 1; i < A.row; i++) {
            solution(i, col) = y(i, col);
            solution(i, col) -= A(i, i - 1) * solution(i - 1, col); //only one subtraction vs. full matrix
            solution(i, col) *= A(i, i);
        }

        //Backsolving. Since tri-diagonal only have to care about one add'l y_i
        for (int i = A.row - 2; i >= 0; i--) {
                solution(i, col) -= A(i, i + 1) * solution(i + 1, col);
        }
    }
    return solution;
}
