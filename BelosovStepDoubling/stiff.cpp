//
// Created by Roberto Fernandez on 11/6/16.
//
// take a simple backward difference step with a single
// Newton iteration
#include "stiff.h"
#include <math.h>
#include <fstream>

void bel(Matrix& fy) {
    Matrix temp(fy);
    fy(0,0) = 77.27 * (temp(1,0) - temp(1,0) * temp(0,0) + temp(0,0) - 8.375e-6 * pow(temp(0,0), 2));
    fy(1,0) = (1 / 77.27) * (-temp(1,0) - temp(1,0)*temp(0,0) + temp(2,0));
    fy(2,0) = 0.161 * (temp(0,0) - temp(2,0));
}

void exp(Matrix& fy) {
    fy(0,0) = -fy(0,0);
}

double RESIDUAL_LIM = 100;
double roundError = 0;
std::ofstream tracef("trace1.0e-2");
std::ofstream results("Belosov Results");
//std::ofstream tracef("ex results");

void step( double time,  double dt,
           Matrix& dy,      // guess at change in y over step
           Matrix& ynew,    // yold+dy
           double& resred,  // residual reduction
           MatrixFP f,      // ode is y' = f
           normFP vnorm ) {  // problem specific norm
    Matrix fa(ynew.row, ynew.col); Matrix Df(ynew.row, ynew.row);
    Matrix ft(f(time, ynew, fa, Df));
    ft = dy - dt * ft; //want ft=0
    Matrix I(ynew.row, ynew.row);
    for(uint i = 0; i < ynew.row; i ++) {
        I(i,i) = 1;
    }
    Df = I - dt * Df;
    factor(Df);
    dy = dy + solve(Df, -ft); //new guess for dy
    ynew = ynew + dy;
    for(uint i = 0; i < ynew.row; i ++) {
        roundError += fa(i,0);
    }

    resred = vnorm(ft)/ (vnorm(dy - dt * f(time, ynew, fa, Df)));
}

// advance from now to the end time, using time
// step control to get a robust solution
void advance( simTime& stime,
              Matrix& Ynow,
              Matrix& Yold,
              MatrixFP f,
              normFP vnorm,
              bool trace ) {
    double ei, resred;
    double maxVal[3] = {INT_MIN, INT_MIN, INT_MIN};
    double maxValt[3] = {0,0,0};
    resred = 0;
    //sets up initial step
    stime.dt = stime.dtOld = stime.dtmin;
    Matrix dmid(Ynow.row, Ynow.col); Matrix dnew(Ynow.row, Ynow.col); Matrix snew(Ynow.row, Ynow.col);
    Matrix dd(Ynow.row, Ynow.col); Matrix ds(Ynow.row, Ynow.col); Matrix s(Ynow.row, Ynow.col);


    while(stime.time < stime.endTime) {
        dd = (stime.dt / (2 * stime.dtOld)) * (Ynow - Yold);
        Matrix ddd(dd);
        dmid = Ynow + dd; //guess for dmid @ tnow + dt/2
        ds = dd + dd;
        snew = dnew = Ynow + ds; //best guess for dnew set to snew also
        step(stime.time, stime.dt, ds, snew, resred, f, vnorm);

        step(stime.time, stime.dt / 2, dd, dmid, resred, f, vnorm);
        step(stime.time, stime.dt / 2, ddd, dnew, resred, f, vnorm);
        ei = vnorm(dnew - snew);

        Matrix YTemp(Ynow + 2 * (dd + ddd) - ds);
        exp(YTemp);
        resred /= vnorm(dd - stime.dt * YTemp);
        roundError *= 1e-15;
        bool residualTest = vnorm(stime.dt * YTemp - ds) < 100 * roundError || resred > RESIDUAL_LIM;
        if ((ei < stime.tol && residualTest) || stime.dt == stime.dtmin) { //accept step
            stime.stepsAccepted++;
            stime.stepsSinceRejection++;
            stime.dtOld = stime.dt;
            stime.time += stime.dt;
            Yold = Ynow;
            Ynow = Ynow + (2 * (dd + ddd) - ds);
            if (ei < 0.25 * stime.tol && stime.stepsSinceRejection > 2)
                stime.dt = std::min(stime.dt * stime.agrow, stime.dtmax);
            else if (ei > 0.75 * stime.tol)
                stime.dt = std::max(stime.dt * stime.ashrink, stime.dtmin);

            if (stime.time + stime.dt > stime.endTime)
                stime.dt = stime.endTime - stime.time;
            else if (stime.time + 2 * stime.dt > stime.endTime)
                stime.dt = (stime.endTime - stime.time) / 2;
            if(stime.time > 500) {
                for(uint i = 0; i < 3; i ++) {
                    if(Ynow(i,0) > maxVal[i]) {
                        maxVal[i] = Ynow(i,0);
                        maxValt[i] = stime.time;
                    }
                }
            }
            std::cout << "time: " << stime.time << "\ty: " << Ynow << std::endl;
            tracef << stime.time << "\t" << Ynow(0, 0) << "\t" << Ynow(1, 0) << "\t" << Ynow(2, 0) << "\n";
        } else {
            stime.stepsRejected++;
            stime.stepsSinceRejection = 0;
            stime.dt = std::max(stime.dt / 2, stime.dtmin);
        }
    }
    results << "time\t value\n";
    for(int i = 0; i < 3; i ++) {
        results << maxValt[i] << '\t' << maxVal[i] << '\n';
    }
    std::cout << "y(" << stime.time << ") = " << Ynow << "Steps Accepted: " <<
              stime.stepsAccepted << "Steps Rejected: " << stime.stepsRejected << std::endl;

}


