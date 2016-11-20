// test full matrix solver

#include "matrix.h"

int main(){

Matrix A(4,4), B(4,2), X(4,2);
unsigned int i,j;
for( i=0; i<4; i++ ){
  for(j=0; j<4; j++)
    A(i,j) = 0.01*i*i+0.02*j;
  A(i,i) += 3+0.3*i;
  for(j=0; j<2; j++)
    X(i,j) = i+1+4*j;
 }
 std::cout << "Test of full matrix\n";
 std::cout << "The results should be:\n"
           << X; //SHOULD BE X

 B = A*X;
 factor( A );
 X = solve(A,B);

std::cout << " And the results are:\n"
          << X;

 std::cout <<"\nTest of tridiagonal matrix\n";
 std::cout << "The results should be:\n"
           << X;

 tridiag C(4,-1.0, 3.0, -1.0 );

 for( i=0; ; i++ ){
   C(i,i) += 0.1*(i+1)-0.001*i*i;
   if( i==3 ) break;
   C(i+1,i) -= 0.02*(i+1);
   C(i,i+1) -=0.03*i*i;
 }

 B = C*X;
 factor(C);
 X = solve(C,B);

std::cout << "The results are:\n"
          << X;

return 0;
}
