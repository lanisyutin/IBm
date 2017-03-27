#include "Calculate_A.h"

void Calculate_A_u(MatrixXd A[5], Grid grid, double Re){

	double d_xx = 1.0 / (grid.d_x*grid.d_x);
	double d_yy = 1.0 / (grid.d_y*grid.d_y);
	int const n1 = grid.N1;
	int const n2 = grid.N2 + 1;


	for (int j = 1; j < (n2 - 1); ++j){
		for (int i = 1; i < (n1 - 1); ++i){
			A[0](i, j) = 1.0 / grid.d_t + (1.0 / Re) * (d_xx + d_yy);
			A[1](i, j) = -1.0 / (2.0*Re*grid.d_y*grid.d_y);
			A[2](i, j) = -1.0 / (2.0*Re*grid.d_x*grid.d_x);
			A[3](i, j) = -1.0 / (2.0*Re*grid.d_y*grid.d_y);
			A[4](i, j) = -1.0 / (2.0*Re*grid.d_x*grid.d_x);

			if (j == 1){
				A[0](i, j) = 1.0 / grid.d_t + (1.0 / Re) * (d_xx + 2.0*d_yy);
				A[1](i, j) = -2.0 / (3.0*Re*grid.d_y*grid.d_y);
				A[2](i, j) = -1.0 / (2.0*Re*grid.d_x*grid.d_x);
				A[3](i, j) = -4.0 / (3.0*Re*grid.d_y*grid.d_y);
				A[4](i, j) = -1.0 / (2.0*Re*grid.d_x*grid.d_x);
			}

			if (j == n2 - 2){
				A[0](i, j) = 1.0 / grid.d_t + (1.0 / Re) * (d_xx + 2.0*d_yy);
				A[1](i, j) = -4.0 / (3.0*Re*grid.d_y*grid.d_y);
				A[2](i, j) = -1.0 / (2.0*Re*grid.d_x*grid.d_x);
				A[3](i, j) = -2.0 / (3.0*Re*grid.d_y*grid.d_y);
				A[4](i, j) = -1.0 / (2.0*Re*grid.d_x*grid.d_x);
			}
		}
	}

	for (int i = 1; i < n1 - 1; ++i){
		int j = 0;
		A[0](i, j) = 1.0;
		A[1](i, j) = 1.0;
		j = n2 - 1;
		A[0](i, j) = 1.0;
		A[1](i, j) = 1.0;


	}

	// outflow du/dx = 0
	for (int j = 0; j < n2; ++j){

		A[0](n1 - 1,j) = 3.0 / (2.0*grid.d_x);
		A[1](n1 - 1,j) = -4.0 / (2.0*grid.d_x);
		A[2](n1 - 1,j) = 1.0 / (2.0*grid.d_x);
		A[0](0,j) = 1.0;
	}



}

void Calculate_A_v(MatrixXd A[5], Grid grid, double Re){

	double d_xx = 1.0 / (grid.d_x*grid.d_x);
	double d_yy = 1.0 / (grid.d_y*grid.d_y);
	int const n1 = grid.N1 + 1;
	int const n2 = grid.N2;


	for (int j = 1; j < (n2 - 1); ++j){
		for (int i = 1; i < (n1 - 1); ++i){
			A[0](i, j) = 1.0 / grid.d_t + (1.0 / Re) * (d_xx + d_yy);
			A[1](i, j) = -1.0 / (2.0*Re*grid.d_y*grid.d_y);
			A[2](i, j) = -1.0 / (2.0*Re*grid.d_x*grid.d_x);
			A[3](i, j) = -1.0 / (2.0*Re*grid.d_y*grid.d_y);
			A[4](i, j) = -1.0 / (2.0*Re*grid.d_x*grid.d_x);

			if (i == 1){
				A[0](i, j) = 1.0 / grid.d_t + (1.0 / Re) * (2.0*d_xx + d_yy);
				A[1](i, j) = -1.0 / (2.0*Re*grid.d_y*grid.d_y);
				A[2](i, j) = -2.0 / (3.0*Re*grid.d_x*grid.d_x);
				A[3](i, j) = -1.0 / (2.0*Re*grid.d_y*grid.d_y);
				A[4](i, j) = -4.0 / (3.0*Re*grid.d_x*grid.d_x);
			}

			if (i == n1 - 2){
				A[0](i, j) = 1.0 / grid.d_t + (1.0 / Re) * (2.0*d_xx + d_yy);
				A[1](i, j) = -1.0 / (2.0*Re*grid.d_y*grid.d_y);
				A[2](i, j) = -4.0 / (3.0*Re*grid.d_x*grid.d_x);
				A[3](i, j) = -1.0 / (2.0*Re*grid.d_y*grid.d_y);
				A[4](i, j) = -2.0 / (3.0*Re*grid.d_x*grid.d_x);
			}

		}
	}

	for (int i = 1; i < n1 - 1; ++i){


		int j = 0;
		A[0](i, j) = -grid.d_y * 0.5;
		A[1](i, j) = 0.0;
		j = n2 - 1;
		A[0](i, j) = grid.d_y * 0.5;
		A[1](i, j) = 0.0;


	}

	// outflow du/dx = 0
	for (int j = 0; j < n2; ++j){

		A[0](n1 - 1,j) = 3.0 / (2.0*grid.d_x);
		A[1](n1 - 1,j) = -4.0 / (2.0*grid.d_x);
		A[2](n1 - 1,j) = 1.0 / (2.0*grid.d_x);
		A[0](0,j) = 1.0;
	}

}

MatrixXd Operator_Ax(MatrixXd A[5], MatrixXd &x, int const n1, int const n2, Grid grid){

	MatrixXd result(n1, n2);
	result.setZero();


	for (int j = 1; j < (n2 - 1); ++j){
		for (int i = 1; i < (n1 - 1); ++i){

			result(i, j) = A[0](i, j) * x(i, j) + A[1](i, j) * x(i, j + 1) + A[2](i, j) * x(i + 1, j) + A[3](i, j) * x(i, j - 1) + A[4](i, j) * x(i - 1, j);
		}
	}

	for (int i = 1; i < n1 - 1; ++i){
		int j = 0;
		result(i, j) = x(i, j);
		j = n2 - 1;
		result(i, j) = x(i, j);
	}

	// outflow du/dx = 0
	for (int j = 0; j < n2; ++j){

		result(n1 - 1, j) = (3.0 * x(n1 - 1, j) - 4.0 * x(n1 - 2, j) + 1.0 * x(n1 - 3, j)) / (2.0*grid.d_x);
		result(0, j) = x(0, j);
	}

	return result;

}