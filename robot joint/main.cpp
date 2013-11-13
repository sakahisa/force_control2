#include <iostream>
#include <cmath>
#include "Eigen/Dense"
#include <vector>

#include "Link.h"
#include "manipulator.h"
#include "constant.h"

using namespace std;
using namespace Eigen;

int main()
{
	Link A;
	manipulator X;
		
	Vector3d x_Lfoot, x_Rfoot;
	Vector3d x_past_Lfoot = Vector3d::Zero(), x_past_Rfoot = Vector3d::Zero();
	Vector3d dx_Lfoot, dx_Rfoot;
	
//	Vector3d x_L2C;
//	Vector3d x_past_L2C = Vector3d::Zero();
//	Vector3d dx_L2C;
	
	Vector3d v_L2R;
	Vector3d v_L2C;
	Vector3d v_COM;
	
//	Vector3d xRef, Error;
	
	VectorXd angles      = VectorXd::Zero(8);
	VectorXd dangles     = VectorXd::Zero(8);
	VectorXd angles_past = VectorXd::Zero(8);
	
	Vector3d dCOMx		= Vector3d::Zero();
	Vector3d COMx		= Vector3d::Zero();
	Vector3d COMx_past	= Vector3d::Zero();
	
//	VectorXd dcomx = Vector3d::Zero();
//	Vector3d error = Vector3d::Zero();
	
	MatrixXd COMJacobian = MatrixXd::Zero(6, angles.size());
	Matrix4d T_L, T_R;
	for(double t = 0.0; t < TMAX + twidth; t += twidth)
	{
		angles << M_PI / 8 * sin(t), -M_PI / 8 * sin(t),
		          M_PI / 8 * sin(t), M_PI / 8 * sin(t),
				  M_PI / 8 * sin(t), M_PI / 8 * sin(t),
		          M_PI / 8 * sin(t), M_PI / 8 * sin(t);
		
		dangles = (angles - angles_past) / twidth;
		angles_past = angles;
		
		T_L = X.forwardKine(-1.0, angles);
		T_R = X.forwardKine(1.0, angles);
		
		x_Lfoot = TRANS(T_L);
		x_Rfoot = TRANS(T_R);
		COMx = (X.COMpos(angles)).head(3);
		
		dx_Lfoot = (x_Lfoot - x_past_Lfoot) / twidth;
		dx_Rfoot = (x_Rfoot - x_past_Rfoot) / twidth;
		dCOMx = (COMx - COMx_past) / twidth;
		
		x_past_Lfoot = x_Lfoot;
		x_past_Rfoot = x_Rfoot;
		COMx_past = COMx;
		
/*		x_L2C = (ROT(T_L)).transpose() * (COMx - x_Lfoot);
		dx_L2C = x_L2C - x_past_L2C;
		x_past_L2C = x_L2C;
*/		
		v_L2R = (ROT(T_L)).transpose() * (dx_Rfoot - dx_Lfoot);
		v_L2C = (ROT(T_L)).transpose() * (dCOMx - dx_Lfoot);
		v_COM = dCOMx;
		
		cout << "test     t = " << t <<endl;
		cout << "test getLfoot2RfootJacobian" << endl;
		X.testJacobian(v_L2R, dangles, X.getLfoot2RfootJacobian(angles));
		cout << "test getLfoot2COMJacobian" << endl;
		X.testJacobian(v_L2C, dangles, X.getLfoot2COMJacobian(angles));
		cout << "test getCOMJacobian" << endl;
		X.testJacobian(v_COM, dangles, X.getCOMJacobian(angles));
		
		}
//	for (int q = 0; q < angles.size() + 1; q ++) 		cout << "X.forwardKine(angles, " << q << ")" << endl << X.forwardKine(angles, q) << endl;
//	for (int p = 0; p < angles.size() + 1; p ++)		cout << X.COMpos(p, angles).transpose() << endl;
	return 0;
}
