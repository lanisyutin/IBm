#pragma once

#include "stdafx.h"
#include "Grid.h"
#include "SolidBody.h"
using namespace Eigen;
using namespace std;

void OutputPressure(MatrixXd data, int n, double output_step, map<int, Circle*> iList, Grid grid);
void OutputVelocity_U(MatrixXd data, int n, int output_step, map<int, Circle*> iList, Grid grid);
void OutputVelocity_V(MatrixXd& data, int n, int output_step, map<int, Circle*> iList, Grid grid);