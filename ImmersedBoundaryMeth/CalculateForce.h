#pragma once

#include "stdafx.h"
#include "Grid.h"
#include "SolidBody.h"

using namespace Eigen;
using namespace std;

double DeltaFunction(double x, double y, Grid grid);
double FunctionD(double r);
void GetInfluenceArea(int& i_min, int& i_max, int& j_min, int& j_max, double x, double y, int size, Grid grid);

double CalculateForce_X(MatrixXd& force_x, map<int, Circle*> &iList, MatrixXd& u, double r, double & Coeff, Grid grid, double alpha_f, double beta_f, double M);
double CalculateForce_Y(MatrixXd& force_y, map<int, Circle*> &iList, MatrixXd& v, double r, double & Coeff, Grid grid, double alpha_f, double beta_f, double M);
