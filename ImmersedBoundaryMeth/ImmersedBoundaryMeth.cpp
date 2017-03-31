#include "stdafx.h"
#include "SolidBody.h"
#include "Grid.h"
#include "Calculate_A.h"
#include "CalculateForce.h"
#include "Calculate_B.h"
#include "Output.h"
#include "BiCGStab.h"
#include "Calculate_press.h"


#pragma warning(disable : 4996)//for using <ctime>
#pragma warning(disable : 4244)//for GetInfluenceArea



using namespace std;
//using namespace Eigen;


// fuctions

void InputData(Grid& grid, double &M, int &Re, double &alpha_f, double &beta_f, double& Zeidel_eps, int& output_step, int& N_max, int& N_Zeidel);
void SetLog(ostream &log, Grid grid, double M, double Re, double alpha_f, double beta_f, double Zeidel_eps);
void PushLog(ostream &log, int n, double eps_u, double eps_v);
void ApplyInitialData(Matrix& u, Grid grid);
void Resize(Matrix& M, const int n, const int m);


int main(){
	int Re;
	int N_max = 0; // number of total iterations
	double alpha_f;
	double beta_f;
	int N_Zeidel; // Number of iterations in Zeidel method
	double Zeidel_eps;
	double Cd; //drag coefficent
	double Cl; // lift coefficent
	double x = 0.0;
	double y = 0.0;
	double r = 0.5;
	double m;
	int output_step = 0; //frequency of output

	// declaration variables
	Grid grid;
	double eps_u = 0.0;
	double eps_v = 0.0;
	double eps_p = 0.0;
	

	int n = 0; // iteration counter
	InputData(grid, m, Re, alpha_f, beta_f, Zeidel_eps,output_step, N_max, N_Zeidel); // Get value of some variables
	CreateMatrix(U_n, grid.N1, grid.N2 + 1);
	CreateMatrix(U_new, grid.N1, grid.N2 + 1);
	CreateMatrix(U_prev, grid.N1, grid.N2 + 1);
	CreateMatrix(B_u, grid.N1, grid.N2 + 1); 
	CreateMatrix(Force_x, grid.N1, grid.N2 + 1); 

	CreateMatrix(V_n, grid.N1+1, grid.N2);
	CreateMatrix(V_new, grid.N1+1, grid.N2);
	CreateMatrix(V_prev, grid.N1 + 1, grid.N2);
	CreateMatrix(B_v, grid.N1 + 1, grid.N2);
	CreateMatrix(Force_y, grid.N1 + 1, grid.N2);
	
	CreateMatrix(P, grid.N1 + 1, grid.N2 + 1);
	CreateMatrix(Delta_P, grid.N1 + 1, grid.N2 + 1);
	CreateMatrix(P_Right, grid.N1 + 1, grid.N2 + 1);
	
	Matrix OperatorA_u[5];
	for (int i = 0; i < 5; i++){
		OperatorA_u[i].resize(grid.N1);
		for (int j = 0; j < grid.N1; j++){
			OperatorA_u[i][j].resize(grid.N2 + 1);
			fill(OperatorA_u[i][j].begin(), OperatorA_u[i][j].end(), 0);
		}
	}

	Matrix OperatorA_v[5];
	for (int i = 0; i < 5; i++){
		OperatorA_v[i].resize(grid.N1+1);
		for (int j = 0; j < grid.N1+1; j++){
			OperatorA_v[i][j].resize(grid.N2);
			fill(OperatorA_v[i][j].begin(), OperatorA_v[i][j].end(), 0);
		}
	}
	

	// list of immersed solids
	map<int, Circle*> solidList;
	ofstream output; // for Drag and Lift coefficents
	ofstream press_output; // press
	ofstream log; // log
	string filename = "Result/coefficent.plt";
	//string filepress = "Result/eps_pressure.plt";
	string filelog = "Result/log.txt";
	log.open(filelog, ios::out);
	SetLog(log,grid,m,Re,alpha_f,beta_f,Zeidel_eps);
	log << endl;

	//нет зависимости начальной скарости от давления
	ApplyInitialData(U_new,grid); // Applying initial data to velocity 
	U_n = U_new;
	U_prev = U_new;

	/*
	Firstly adding some circles
	*/
	Circle c1(3.5, 2.1, r, n, grid);
	Circle c2(3.5, 4.9, r, n, grid);
	Circle c3(1.5, 1.9, r, n, grid);
	Circle c4(1.5, 5.1, r, n, grid);
	c1.AddSolid(solidList);
	c2.AddSolid(solidList);
	c3.AddSolid(solidList);
	c4.AddSolid(solidList);
	


	CalculateForce_X(Force_x, solidList, U_new, r, Cd, grid, alpha_f, beta_f,m);
	CalculateForce_Y(Force_y, solidList, V_new, r, Cl, grid, alpha_f, beta_f,m);

	OutputVelocity_U(U_new, -1, output_step, solidList, grid);
	OutputVelocity_V(V_new, -1,	output_step, solidList, grid);
	//press_output.open(filepress);
	while (n <= N_max){
		//creation new solids

		if (n > 0 && fmod((double)n*grid.d_t, 1.5) == 0.0){
			int chance = 80;
			if (rand() % 100 + 1 <= chance){
				x =1 + ((rand() % 100 + 1) / 100.0);//ПРОВЕРИТЬ ЧТО РЕЗУЛЬТАТ DOUBLE
				y =1 + ((rand() % 200 + 1) / 100.0);
				Circle c(x,y, r, n, grid);
				c.AddSolid(solidList);
			}

			if (rand() % 100 + 1 <= chance){
				x =1 + ((rand() % 100 + 1) / 100.0);
				y =4 + ((rand() % 200 + 1) / 100.0);
				Circle c(x, y, r, n, grid);
				c.AddSolid(solidList);
			}

			if (rand() % 100 + 1 <= chance){
				x =3+ ((rand() % 100 + 1) / 100.0);
				y =1+ ((rand() % 200 + 1) / 100.0);
				Circle c(x, y, r, n, grid);
				c.AddSolid(solidList);
			}

			if (rand() % 100 + 1 <= chance){
				x =3+ ((rand() % 100 + 1) / 100.0);
				y =4+ ((rand() % 200 + 1) / 100.0);
				Circle c(x, y, r, n, grid);
				c.AddSolid(solidList);
			}
		}

		eps_u = 0.0;
		eps_v = 0.0;

		B_u = CalculateB_u(U_n, V_n, U_prev, V_prev, P, Force_x, grid, Re);
		B_v = CalculateB_v(U_n, V_n, U_prev, V_prev, P, Force_y, grid, Re);

		Calculate_A_u(OperatorA_u, grid, Re);//
		Calculate_A_v(OperatorA_v, grid, Re);//

		BiCGStab(U_new, grid.N1, grid.N2 + 1, OperatorA_u, B_u,grid);
		BiCGStab(V_new, grid.N1 + 1, grid.N2, OperatorA_v, B_v,grid);

		P_Right = Calculate_Press_Right( U_n, V_n,grid);

		for (int i = 0; i < Delta_P.size(); ++i)
			for (int j = 0; j < Delta_P[i].size(); ++j){
			Delta_P[i][j] = 0.0;
			}
		

		eps_p = Calculate_Press_correction(Delta_P, P_Right,N_Zeidel,Zeidel_eps,grid); //Должна ли поправка зависеть от скорости?
		press_output << n << ' ' << eps_p << endl;



		for (int i = 0; i < grid.N1 + 1; ++i){
			for (int j = 0; j < grid.N2 + 1; ++j){
				P[i][j] = P[i][j] + 0.8 * Delta_P[i][j];
			}
		}

		for (int i = 1; i < grid.N1 - 1; ++i){
			for (int j = 1; j < grid.N2; ++j){
				U_new[i][j] = U_new[i][j] - grid.d_t * (Delta_P[i + 1][j] - Delta_P[i][j]) / grid.d_x;
			}
		}

		for (int j = 1; j < grid.N2; ++j){
			int i = grid.N1 - 1;
			U_new[i][j] = U_new[i - 1][j];
		}


		for (int i = 1; i < grid.N1 + 1; ++i){
			for (int j = 1; j < grid.N2 - 1; ++j){
				V_new[i][j] = V_new[i][j] - grid.d_t * (Delta_P[i][j+1] - Delta_P[i][j]) / grid.d_y;
			}
		}
		
		for (int i = 0; i < grid.N1; ++i){
			for (int j = 0; j < grid.N2 + 1; ++j){
				if (fabs(U_n[i][j] - U_new[i][j]) > eps_u){
					eps_u = fabs(U_n[i][j] - U_new[i][j]);
				}

				U_prev[i][j] = U_n[i][j];
				U_n[i][j] = U_new[i][j];
			}
		}

		for (int i = 0; i < grid.N1 + 1; ++i){
			for (int j = 0; j < grid.N2; ++j){
				if (fabs(V_n[i][j] - V_new[i][j]) > eps_v){
					eps_v = fabs(V_n[i][j] - V_new[i][j]);
				}
				V_prev[i][j] = V_n[i][j];
				V_n[i][j] = V_new[i][j];
			}
		}


		CalculateForce_X(Force_x, solidList, U_new, r, Cd, grid,alpha_f,beta_f,m);
		CalculateForce_Y(Force_y, solidList, V_new, r, Cl, grid, alpha_f, beta_f,m);

		
		for (auto& solid = solidList.begin(); solid != solidList.end(); ++solid){
			
			if ((solid->second->moveSolid == false)&&(n - solid->second->start_n > 10)){
				solid->second->moveSolid = true;
			}

			if (solid->second->moveSolid){
				//update position
				for (int k = 0; k < grid.NF; ++k){
					solid->second->Bound[0][k] += solid->second->U * grid.d_t;
					solid->second->Bound[1][k] += solid->second->V * grid.d_t;
					solid->second->x += solid->second->U * grid.d_t;
					solid->second->y += solid->second->V * grid.d_t;
				}
			}
			//delete bodies which move 80% of length
			if (solid->second->x > grid.L*0.8){
					solidList.erase(solid->first);
				
				}

		}
		//

		if (n < 1000 || (n > 1000 && 0 == n % 100)){
			cout << "n  = " << n << " | eps_u = " << eps_u << " | eps_v = " << eps_v << "		";
			time_t t = chrono::system_clock::to_time_t(chrono::system_clock::now());   // get time now
			cout << ctime(&t);

			PushLog(log, n, eps_u, eps_v);
			log.flush();
		}

		if (0 == n % output_step){
			OutputVelocity_U(U_new, n, output_step, solidList, grid);
			OutputVelocity_V(V_new, n, output_step, solidList, grid);
			OutputPressure(P, n,output_step, solidList,grid);
		}


#pragma message("Clean hardcode for tolerance")
		if (eps_u < 1e-3 && eps_v < 1e-3){
			OutputVelocity_U(U_new, n, output_step, solidList, grid);
			OutputVelocity_V(V_new, n, output_step, solidList, grid);
			OutputPressure(P, n, output_step, solidList, grid);
			break;
		}


		++n;
	}
	log.close();
	cout << "Over" << endl;
	getchar();

	return 0;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void SetLog(ostream& log, Grid grid, double M, double Re, double alpha_f, double beta_f, double Zeidel_eps){

	log << "The IBM program starts.		";
	time_t t = chrono::system_clock::to_time_t(chrono::system_clock::now());   // get time now
	log << ctime(&t) << endl;
	log << "The parameters are the following:"<<endl;
	log << "Mass of a particle            : M   = " << M << endl;
	log << "Reynolds number               : Re  = " << Re << endl;
	log << "Channel length                : L   = " << grid.L << endl;
	log << "Channel width                 : W   = " << grid.H << endl;
	log << "Number of nodes on            : Nх  = " << grid.N1 << endl;
	log << "Number of nodes on            : Ny  = " << grid.N2<< endl;
	log << "Number of nodes for a particle: Np  = " << grid.NF<< endl;
	log << "Time step                     : tau = " << grid.d_t<< endl;
	log << "Force parameter alpha         : alpha = " <<alpha_f<< endl;
	log << "Force parameter beta          : beta  = " <<beta_f << endl;
	log << "Tolerance for Zeidel method   : tol = " << Zeidel_eps << endl;

}

void PushLog(ostream& log, int n, double eps_u, double eps_v){
	time_t t = chrono::system_clock::to_time_t(chrono::system_clock::now()); 
	log << "n  = " << n << " | eps_u = " << eps_u << " | eps_v = " << eps_v << '\t';
	log << ctime(&t);
}

// Apply initial data for velocity
void ApplyInitialData(Matrix &u, Grid grid){

	// Poiseuille flow 
	for (int i = 0; i < grid.N1; ++i){
		for (int j = 1; j < grid.N2; ++j){
			u[i][j] = (pow((grid.H) / 2.0, 2) - pow((j - 0.5)*grid.d_y - grid.H/ 2.0, 2));
		}
	}
}

void InputData(Grid& grid, double &M, int &Re, double &alpha_f, double &beta_f, double& Zeidel_eps, int& output_step, int& N_max, int& N_Zeidel){

	ifstream input;
	string filename = "input.txt";
	input.open(filename.c_str());

	input >> M;
	input >> Re;
	input >> grid.L;
	input >> grid.H;
	input >> grid.N1;
	input >> grid.N2;
	input >> grid.NF;
	input >> grid.d_t;
	input >> alpha_f;
	input >> beta_f;
	input >> output_step;
	input >> N_max;
	input >> N_Zeidel;
	input >> Zeidel_eps;
	input.close();

	grid.d_x = grid.L / (grid.N1 - 1);
	grid.d_y = grid.H / (grid.N2 - 1);

}

void Resize(Matrix& M, const int n, const int m){
	for (int i = 0; i < n; i++){
		
	}
}