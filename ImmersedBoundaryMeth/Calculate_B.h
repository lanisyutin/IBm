#pragma once
#include "stdafx.h"
#include "Grid.h"
using namespace Eigen;
using namespace std;


MatrixXd CalculateB_u(MatrixXd &u_n, MatrixXd &v_n, MatrixXd &u_prev, MatrixXd &v_prev, MatrixXd &p, MatrixXd &force, Grid grid, double Re);
MatrixXd CalculateB_v(MatrixXd &u_n, MatrixXd &v_n, MatrixXd &u_prev, MatrixXd &v_prev, MatrixXd &p, MatrixXd &force, Grid grid, double Re);
