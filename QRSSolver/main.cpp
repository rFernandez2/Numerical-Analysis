#include <iostream>
#include "qrsolve.h"

int main() {
    Matrix A(4,4), B(4,2), X(4,2), Q(4,4), R(4,4);
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
              << A; //SHOULD BE X


    factor( A, Q, R );

    std::cout << " And the results are:\n"
              << A << '=' << Q << "*" << R << std::endl;
    Matrix ans(Q * R);
    std::cout << "so we have our answer:\n" << ans << std::endl;
    return 0;
}