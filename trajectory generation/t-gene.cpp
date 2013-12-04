#include <iostream>
#include <cmath>
#include <fstream>
#include "Eigen/Dense"
#include <vector>

#define twidth 0.005
#define N 900
#define g 9.8
#define Hcom 0.5

using namespace std;
using namespace Eigen;


int main()
{
	ofstream ofs( "t-gene.txt" );
	Matrix3d A;
	A(0, 0) = 1;		A(0, 1) = twidth;		A(0, 2) = (twidth * twidth) / 2;
	A(1, 0) = 0;		A(1, 1) = 1;		A(1, 2) = twidth;
	A(2, 0) = 0;		A(2, 1) = 0;		A(2, 2) = 1;
	
	Vector3d B;
	B << (twidth * twidth * twidth) / 6, (twidth * twidth) / 2, twidth;
	
	Vector3d C;
	C << 1, 0, -Hcom / g;
	
	MatrixXd PX = MatrixXd::Zero(N ,3);
	Vector3d x0 = Vector3d::Zero();
	VectorXd zmpRef = VectorXd::Zero(N);
	
	for (int i = N / 2; i < N; i++)	zmpRef(i) = 0.1;
	
	for (int j = 0; j < N; j++)
	{
		if(j ==0)
			PX.row(j) = C.transpose() * A;
		else
			PX.row(j) = PX.row(j - 1) * A;
	}
	
	MatrixXd PU = MatrixXd::Zero(N, N);
	Matrix3d AA;
	for (int k = 0; k < N; k++)
	{
		AA = Matrix3d::Identity();
		
		for (int l = k; l > 0; l--)
		{
			PU(k, l) = C.transpose() * AA * B;
			AA *= A;
		}
	}
	
	double Q = 1.0;
	double R = 1E-6;
	MatrixXd I = MatrixXd::Identity(N, N);
	MatrixXd w = Q * PU.transpose() * PU + R * I;
	VectorXd v = Q * PU.transpose() * PX * x0 - Q * PU.transpose() * zmpRef;
	VectorXd jerk = -w.inverse() * v;
	
	Vector3d xk = x0;
	double zk = 0.0;
	double t = 0.0;
	for (int m = 1; m < N; m++)
	{
		xk = A * xk + B * jerk(m);
		zk = C.transpose() * xk;
		t = twidth * m;
		
		ofs << t << "	" << xk(0) << "	" << zk << endl;
	}
}
