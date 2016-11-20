// test of 3rd-order Adams P-C

#include "utility.h"
#include "numvec.h"
#include "adams.h"

#include <cmath>
#include <ostream>

NumVec fpoly( double time, const NumVec& Y );
NumVec fcomet( double time, const NumVec& Y );


NumVec fpoly( double time, const NumVec& Y )
{// time^j, j=0,1,2,3;
  NumVec val(4);
  val[0] = 1;
  val[1] = 2 * time;
  val[2] = 3 * time * time;
  val[3] = 4 * time * time * time;
  return val;
}

NumVec fcomet( double time, const NumVec& Y )
{
  NumVec val(4);
  val[0] = Y[2];
  val[1] = Y[3];
  double rs = sq(Y[0]) +sq(Y[1]); // r*r
  double rc = rs*sqrt(rs);        // r*r*r
  val[2] = -Y[0]/rc;
  val[3] = -Y[1]/rc;
  return val;

}

void testComet(double tol) {
  simTime stime;
  double s = 0.3;
  NumVec Ynow(4), Yold(4);
  Ynow[0] = 1.0;
  Ynow[1] = 0.0;
  Ynow[2] = 0.0;
  Ynow[3] = s;
  double period = 2*M_PI/pow(2.-s*s,1.5);
  stime.endTime = 3*period;
  stime.tol = tol;
  stime.agrow = 1.25;
  stime.ashrink = 0.8;
  stime.stepsSinceRejection = 0;
  stime.stepsAccepted = 0;
  stime.stepsRejected = 0;
  double
          energy = 0.5*(sq(Ynow[2])+sq(Ynow[3]))
                   - 1.0/sqrt(sq(Ynow[0])+sq(Ynow[1]));
  std::cout << "At the start time " << stime.time
            << " the solution is\n"  << Ynow
            << "The energy is " << energy
            << ".\n";

  advance( stime, Ynow, Yold, fcomet, true);

  energy = 0.5*(sq(Ynow[2])+sq(Ynow[3]))
           - 1.0/sqrt(sq(Ynow[0])+sq(Ynow[1]));
  std::cout << "At the end time " << stime.time
            << " the solution is\n" << Ynow
            << "The energy is " << energy
            << ".\n";
  std::cout << "The number of accepted steps is "
            << stime.stepsAccepted << ".\n";
  std::cout << "The number of rejected steps is "
            << stime.stepsRejected <<".\n";
  std::cout << "The tolerance is " << stime.tol << ".\n";
}

void testPoly(double tol) {
  simTime stime;
  NumVec Ynow(4), Yold(4);
  Ynow[0] = Ynow[1] = Ynow[2] = Ynow[3] = 0;
  stime.endTime = 2;
  stime.tol = tol;
  stime.agrow = 1.25;
  stime.ashrink = 0.8;
  stime.stepsSinceRejection = 0;
  stime.stepsAccepted = 0;
  stime.stepsRejected = 0;
  std::cout << "At the start time " << stime.time
            << " the solution is\n"  << Ynow;

  advance( stime, Ynow, Yold, fpoly, false);

  std::cout << "At the end time " << stime.time
            << " the solution is\n" << Ynow;
  std::cout << "The number of accepted steps is "
            << stime.stepsAccepted << ".\n";
  std::cout << "The number of rejected steps is "
            << stime.stepsRejected <<".\n";
  std::cout << "The tolerance is " << stime.tol << ".\n";
}

int main()
{
  testComet(0.01);
  testComet(0.001);
  testPoly(0.01);
  return 0;
}



