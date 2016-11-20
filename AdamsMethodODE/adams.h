// interface to a 3rd order Adams-Bashforth solver
#include "numvec.h"

typedef NumVec (*NumVecFP)( double , const NumVec&  );

class simTime{ // class for time step control
 public:
  double time;
  double dt;
  double dtOld;
  double tol;
  double agrow;
  double ashrink;
  double dtmin;
  double dtmax;
  double endTime;
  int stepsSinceRejection;
  int stepsRejected;
  int stepsAccepted;
  simTime()
    :time(0.0), dt(0.0), dtOld(0.0), tol(1.0e-3), agrow(1.25),
    ashrink(1.0/1.25), dtmin(1.0e-6), dtmax(1.0), endTime(1.0), 
    stepsSinceRejection(0), stepsRejected(0), stepsAccepted(0){}
};

void step( double time,  double dt, double dtold,
           const NumVec& fold,
           const NumVec& fnow,
           const NumVec& Ynow,
           NumVec& Ynew,
           NumVec& ddy,
           NumVecFP f );

void advance( simTime& stime, NumVec& Ynow, NumVec& Yold, NumVecFP f, bool trace=true );
