#include <stdio.h>
#include <string>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstring>
#include "numvec.h"
#include "utility.h"
#include "adams.h"

std::ofstream tracef("trace1.0e-2");
std::ofstream traceff("trace1.0e-3");
std::ofstream tracepoly("fpoly");

void step(double time,  double dt, double dtold, const NumVec& fold, const NumVec& fnow,
          const NumVec& Ynow, NumVec& Ynew, NumVec& ddy, NumVecFP f ) {
    double tnew = time + dt;
    NumVec l = (1/ dtold) * (fnow - fold);
    //Ynew = Ynow + (1 / dtold) * (fold + fnow) + (1/dt) * (f(time, Ynew) + fnow);
    //std::cout << "Ynow: " << Ynow << "fnow: " << fnow << std::endl;
    NumVec lTnew = fnow + (dt / dtold) * (fnow - fold); //linear extrapolation of approx now value
    NumVec w = Ynow + (dt * fnow) + (dt * dt / 2) * l; //prediction
    NumVec fguess = f(tnew, w);
    double num = dt * (2 * (tnew) + time - 3 * (time - dtold));
    double den = 6 * (dt);
    Ynew = w + (num / den) *  (fnow - lTnew);
    ddy = Ynew - w;
}

//Constant step-size function that works
/*void advance(simTime& stime, NumVec& Ynow, NumVec& Yold, NumVecFP f, bool trace) {
    stime.dt = (stime.endTime - stime.time) / 10000;
    stime.dtOld = stime.dt;
    NumVec fold = f(stime.time, Ynow);
    NumVec err;
    for(int i = 0; i < Yold.size(); i ++) { //there has to be a better way to do this
        Yold[i] = Ynow[i];
    }
    while(stime.time < stime.endTime) {
        NumVec fnow = f(stime.time, Ynow);
        step(stime.time, stime.dt, stime.dtOld, fold, fnow, Yold,
             Ynow, err, f);
        fold = fnow;
        for(int i = 0; i < Yold.size(); i ++) { //there has to be a better way to do this
            Yold[i] = Ynow[i];
        }
        stime.time += stime.dt;
    }
}*/

void advance(simTime& stime, NumVec& Ynow, NumVec& Yold, NumVecFP f, bool trace) {
    stime.dt = stime.dtOld = stime.dtmin;

    NumVec err;
    NumVec fold = f(stime.time, Ynow);
    double ei;

    Yold = Ynow;

    while(stime.time <= stime.endTime) {
        NumVec fnow = f(stime.time, Ynow);
        step(stime.time, stime.dt, stime.dtOld, fold, fnow, Yold,
             Ynow, err, f);
        ei = norm(err); //calculates euclidean norm on error matrix
        if(ei > stime.tol && 2 * stime.dt > stime.dtmin) { //reject step and rollback
            stime.stepsSinceRejection = 0;
            stime.stepsRejected ++;
            stime.dt /= 2;
            continue;
        } else { //accept step
            if(trace && stime.tol == 0.01) {
                //write comet results
                tracef << stime.time << "\t" <<  Ynow << "\n";
            } else if(trace && stime.tol == 0.001) {
                traceff << stime.time << "\t" <<  Ynow[2] << "\n";
            } else {
                tracepoly << stime.time << "\t" << Ynow << "\n";
            }
            Yold = Ynow;
            fold = fnow;
            stime.stepsSinceRejection ++;
            stime.stepsAccepted ++;
            stime.dtOld = stime.dt;
            stime.time += stime.dt;
            if(stime.stepsSinceRejection > 1 && ei < 4 * stime.tol)
                stime.dt *= stime.agrow;
            else if(ei > 0.75 * stime.tol)
                stime.dt *= stime.ashrink;
            if(stime.time + stime.dt > stime.endTime)
                stime.dt = stime.endTime - stime.time;
            else if(stime.time + 2 * stime.dt > stime.endTime)
                stime.dt = (stime.endTime - stime.time) / 2;
            if(!(stime.time < stime.endTime && stime.endTime - stime.time < 2 * stime.dtmin)) {
                if(stime.dt > stime.dtmax)
                    stime.dt = stime.dtmax;
                else if(stime.dt < stime.dtmin)
                    stime.dt = stime.dtmin;
            }
        }

        //std::cout << "Time: " << stime.time << "\t Step: " << stime.dt
        //<< "\tdtold: " << stime.dtOld << "\tYnow: " << Ynow << std::endl;
    }
}

