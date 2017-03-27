#pragma once

#include "stdafx.h"
#include "Grid.h"
#include "Calculate_A.h"
using namespace Eigen;
using namespace std;


/* We want to solve Poisson equation for velocity.
Au = b, where A MatrixXd of coefficents of leap-frog scheme applyied to poisson equation. and b right side
preparation for CG
Suppose in u first approximation ( in fact in u - velocity fromprevious step)
r(0) = b - Au
z(0) = r(0)
in b_norm calculate Euclid norm of vector b*/
void BiCGStab(MatrixXd& res, int const n1, int const n2, MatrixXd operator_A[5], MatrixXd &b,Grid grid);
double ScalarOperator(MatrixXd &a, MatrixXd &b, int const n1, int const n2);
