#pragma once

#include "stdafx.h"
#include "Grid.h"
using namespace Eigen;
using namespace std;






void Calculate_A_u(MatrixXd A[5], Grid grid, double Re);
void Calculate_A_v(MatrixXd A[5], Grid grid, double Re);
MatrixXd Operator_Ax(MatrixXd A[5], MatrixXd &x, int n1, int n2, Grid grid);