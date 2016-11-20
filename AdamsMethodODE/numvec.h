#ifndef NUMVECH
#define NUMVECH
// interface for numerical vectors

#include <cstdlib>
#include <iostream>
#include<vector>

#include "utility.h"

typedef std::vector<double> NumVec;


NumVec operator+(const NumVec& A, const NumVec& B); // A+B
NumVec operator*( double a, const NumVec& B); // a*B
NumVec operator-(const NumVec& A, const NumVec& B); // A-B
double dot(const NumVec &A, const NumVec &B); //A,B
double norm(const NumVec &A);

std::ostream& operator<<(std::ostream& os, const NumVec& A);
#endif






