#define _USE_MATH_DEFINES
#include <math.h>
#include "CalculateForce.h"

double DeltaFunction(double x, double y, Grid grid){
	return 1.0 / (grid.d_x*grid.d_y) * FunctionD(x / grid.d_x) * FunctionD(y / grid.d_y);
}

double FunctionD(double r){
	if ((0.0 <= fabs(r)) && (fabs(r) < 1.0)){
		return 1.0 / 8.0*(3.0 - 2.0 * fabs(r) + sqrt(1.0 + 4.0 * fabs(r) - 4.0 * r * r));
	}
	if ((1.0 <= fabs(r)) && (fabs(r) < 2.0)){
		return 1.0 / 8.0*(5.0 - 2.0 * fabs(r) - sqrt(-7.0 + 12.0 * fabs(r) - 4.0 * r * r));
	}
	if (2.0 <= fabs(r)){
		return 0.0;
	}
	return 0;
}
void GetInfluenceArea(int& i_min, int& i_max, int& j_min, int& j_max, double x, double y, int size, Grid grid){

	i_max = (x / grid.d_x) + size;
	i_min = (x / grid.d_x) - size;

	j_max = (y / grid.d_y) + size;
	j_min = (y / grid.d_y) - size;

	if (i_min < 0){
		i_min = 0;
	}
	if (j_min < 0){
		j_min = 0;
	}
}

double CalculateForce_X(MatrixXd& force_x, map<int, Circle*> &iList, MatrixXd& u, double r, double & Coeff, Grid grid, double alpha_f, double beta_f, double M){

	int const n1 = grid.N1, const n2 = grid.N2+1;

	vector<double> bound_Force_x;
	double bound_norm = 0.0;
	double new_x = 0.0;
	double new_u_bound = 0.0;
	double f1 = 0.0;

	bound_Force_x.resize(grid.NF);

	/*
	calculating force F for Lagrange
	in new_x and new_y calculated value of velocity in Lagrangian point. It calculates by using near points and discrete delta function
	*/

	for (auto& solid : iList){
		MatrixXd force_x_temp(n1, n2);
		force_x_temp.setZero();
		for (int k = 0; k < grid.NF; ++k){

			new_x = 0.0;
			int i_max = 0;
			int i_min = 0;
			int j_max = 0;
			int j_min = 0;

			GetInfluenceArea(i_min, i_max, j_min, j_max, solid.second->Bound[0][k], solid.second->Bound[1][k], 3, grid);

			if (i_max >= n1){
				i_max = n1 - 1;
			}
			if (j_max >= n2){
				j_max = n2 - 1;
			}

			for (int i = i_min; i <= i_max; ++i){
				for (int j = j_min; j <= j_max; ++j){

					new_x += u(i, j) * DeltaFunction(i*grid.d_x - solid.second->Bound[0][k], (j - 0.5)*grid.d_y - solid.second->Bound[1][k],grid) * grid.d_x * grid.d_y;

				}
			}


			solid.second->Integral_x[k] += (new_x - solid.second->U) * grid.d_t;
			bound_Force_x[k] = alpha_f * solid.second->Integral_x[k] + beta_f * (new_x - solid.second->U);
		}

		//calculating force f for Euler points
		// spreading force by delta function
		for (int k = 0; k < grid.NF; ++k){

			int i_max = 0;
			int i_min = 0;

			int j_max = 0;
			int j_min = 0;

			GetInfluenceArea(i_min, i_max, j_min, j_max, solid.second->Bound[0][k], solid.second->Bound[1][k], 3,grid);

			if (i_max >= n1){
				i_max = n1 - 1;
			}
			if (j_max >= n2){
				j_max = n2 - 1;
			}

			for (int i = i_min; i <= i_max; ++i){
				for (int j = j_min; j <= j_max; ++j){

					force_x_temp(i, j) += bound_Force_x[k] * DeltaFunction(i*grid.d_x - solid.second->Bound[0][k], (j - 0.5)*grid.d_y - solid.second->Bound[1][k], grid) * solid.second->d_s * solid.second->d_s;
				}
			}
		}


		int i_max = 0;
		int i_min = n1;

		int j_max = 0;
		int j_min = n2;

		for (int k = 0; k < grid.NF; ++k){

			int i_max_temp = 0;
			int i_min_temp = 0;

			int j_max_temp = 0;
			int j_min_temp = 0;

			GetInfluenceArea(i_min_temp, i_max_temp, j_min_temp, j_max_temp, solid.second->Bound[0][k], solid.second->Bound[1][k], 3, grid);

			if (i_max_temp > i_max){
				i_max = i_max_temp;
			}
			if (i_min_temp < i_min){
				i_min = i_min_temp;
			}
			if (j_max_temp > j_max){
				j_max = j_max_temp;
			}
			if (j_min_temp < j_min){
				j_min = j_min_temp;
			}

		}

		if (i_max >= n1){
			i_max = n1 - 1;
		}
		if (j_max >= n2){
			j_max = n2 - 1;
		}

		Coeff = 0.0;
		for (int i = i_min; i <= i_max; ++i){
			for (int j = j_min; j <= j_max; ++j){


				force_x(i, j) += force_x_temp(i, j);

				Coeff += force_x_temp(i, j) * grid.d_x * grid.d_y;

			}

		}

		if (solid.second->moveSolid){

			solid.second->U = solid.second->U + (-Coeff * grid.d_t) / (M - M_PI * solid.second-> r * solid.second->r);
		}
	}

	return 0;

}

double CalculateForce_Y(MatrixXd& force_y, map<int, Circle*> &iList, MatrixXd& v, double r, double & Coeff, Grid grid, double alpha_f, double beta_f, double M){

	int const n1 = grid.N1+1, const n2 = grid.N2 ;
	vector<double> bound_Force_y;

	double bound_norm = 0.0;

	double new_y = 0.0;

	double new_v_bound = 0.0;

	MatrixXd force_y_temp(n1,n2);
	force_y_temp.setZero();
	bound_Force_y.resize(grid.NF);


	/*
	for (int i = 0; i < n1; ++i){
		for (int j = 0; j < n2; ++j){
			force_y(i, j) = 0.0;
		}
	}
	*/

	/*calculating force F for Lagrange

	in new_x and new_y calculated value of velocity in Lagrangian point. It calculates by using near points and discrete delta function

	*/

	for (auto& solid : iList){
		/*
		for (int i = 0; i < n1; ++i){
			for (int j = 0; j < n2; ++j){
				force_y_temp(i, j) = 0.0;
			}
		}
		*/
		for (int k = 0; k < grid.NF; ++k){

			new_y = 0.0;

			int i_max = 0;
			int i_min = 0;

			int j_max = 0;
			int j_min = 0;

			GetInfluenceArea(i_min, i_max, j_min, j_max, solid.second->Bound[0][k], solid.second->Bound[1][k], 3,grid);

			if (i_max >= n1){
				i_max = n1 - 1;
			}
			if (j_max >= n2){
				j_max = n2 - 1;
			}

			for (int i = i_min; i <= i_max; ++i){
				for (int j = j_min; j <= j_max; ++j){


					new_y += v(i, j) * DeltaFunction((i - 0.5)*grid.d_x - solid.second->Bound[0][k], j*grid.d_y - solid.second->Bound[1][k],grid) * grid.d_x * grid.d_y;

				}
			}


			solid.second->Integral_y[k] += (new_y - solid.second->V) * grid.d_t;

			bound_Force_y[k] = alpha_f * solid.second->Integral_y[k] + beta_f * (new_y - solid.second->V);

		}



		//calculating force f for Euler points
		// spreading force by delta function
		for (int k = 0; k < grid.NF; ++k){

			int i_max = 0;
			int i_min = 0;

			int j_max = 0;
			int j_min = 0;

			GetInfluenceArea(i_min, i_max, j_min, j_max, solid.second->Bound[0][k], solid.second->Bound[1][k], 3,grid);

			if (i_max >= n1){
				i_max = n1 - 1;
			}
			if (j_max >= n2){
				j_max = n2 - 1;
			}

			for (int i = i_min; i <= i_max; ++i){
				for (int j = j_min; j <= j_max; ++j){

					force_y_temp(i, j) += bound_Force_y[k] * DeltaFunction((i - 0.5)*grid.d_x - solid.second->Bound[0][k], j*grid.d_y - solid.second->Bound[1][k], grid) * solid.second->d_s * solid.second->d_s;
				}
			}

		}

		int i_max = 0;
		int i_min = n1;

		int j_max = 0;
		int j_min = n2;

		for (int k = 0; k < grid.NF; ++k){

			int i_max_temp = 0;
			int i_min_temp = 0;

			int j_max_temp = 0;
			int j_min_temp = 0;

			GetInfluenceArea(i_min_temp, i_max_temp, j_min_temp, j_max_temp, solid.second->Bound[0][k], solid.second->Bound[1][k], 3,grid);

			if (i_max_temp > i_max){
				i_max = i_max_temp;
			}
			if (i_min_temp < i_min){
				i_min = i_min_temp;
			}
			if (j_max_temp > j_max){
				j_max = j_max_temp;
			}
			if (j_min_temp < j_min){
				j_min = j_min_temp;
			}
		}

		if (i_max >= n1){
			i_max = n1 - 1;
		}
		if (j_max >= n2){
			j_max = n2 - 1;
		}

		Coeff = 0.0;
		for (int i = i_min; i <= i_max; ++i){
			for (int j = j_min; j <= j_max; ++j){


				force_y(i, j) += force_y_temp(i, j);

				Coeff += force_y_temp(i, j) * grid.d_x * grid.d_y;

			}

		}

		if (solid.second->moveSolid){

			solid.second->V = solid.second->V + (-Coeff * grid.d_t) / (M - M_PI * solid.second->r * solid.second->r);
		}
	}


	return 0;

}