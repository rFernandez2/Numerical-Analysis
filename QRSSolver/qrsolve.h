//
// Created by Roberto Fernandez on 11/16/16.
//
#include "utility.h"
#include "matrix.h"
#include <math.h>

#ifndef QRSSOLVER_QRSOLVE_H
#define QRSSOLVER_QRSOLVE_H

void factor(Matrix &A, Matrix& Q, Matrix& R); //factors using householder relfections

Matrix solve(Matrix &A, Matrix &y); //solves assuming factored with HH reflections

#endif //QRSSOLVER_QRSOLVE_H
