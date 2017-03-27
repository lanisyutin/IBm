#pragma once

#include "stdafx.h"
#include "Grid.h"
using namespace Eigen;
using namespace std;


double Calculate_Press_correction(MatrixXd& delta_p, MatrixXd &b_p, MatrixXd &u, int const N_Zeidel, double const Zeidel_eps, Grid grid);
MatrixXd Calculate_Press_Right(MatrixXd& u, MatrixXd& v, Grid grid);