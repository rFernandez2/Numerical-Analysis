// interface to a 2nd order backward difference method
#include "matrix.h"

typedef Matrix (*MatrixFP)( double t, const Matrix& y, Matrix& fa, Matrix& df  );
typedef double (*normFP)( const Matrix& y );

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
            :time(0.0), dt(0.0), dtOld(0.0), tol(1.0e-3),
             agrow(2), ashrink(1.0/1.25),
             dtmin(1.0e-3),dtmax(10.0), endTime(1.0),
             stepsSinceRejection(0), stepsRejected(0), stepsAccepted(0){}
};

// take a simple backward difference step with a single
// Newton iteration
void step( double time,  double dt,
           Matrix& dy,      // guess at change in y over step
           Matrix& ynew,    // yold+dy
           double& resred,  // residual reduction
           MatrixFP f,      // ode is y' = f
           normFP vnorm,
           double& roundError );  // problem specific norm

// advance from now to the end time, using time
// step control to get a robust solution
void advance( simTime& stime,
              Matrix& Ynow,
              Matrix& Yold,
              MatrixFP f,
              normFP vnorm,
              bool trace=true );
