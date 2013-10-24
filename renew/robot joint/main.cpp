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

//	Vector3d xRef, Error;
	
	VectorXd angles      = VectorXd::Zero(8);
	VectorXd dangles     = VectorXd::Zero(8);
	VectorXd angles_past = VectorXd::Zero(8);
	
//	Vector4d dCOMx		= Vector4d::Zero();
	Vector4d COMx		= Vector4d::Zero();
//	Vector4d COMx_past	= Vector4d::Zero();
	
//	VectorXd dcomx = Vector3d::Zero();
//	Vector3d error = Vector3d::Zero();
	
	MatrixXd COMJacobian = MatrixXd::Zero(6, angles.size());
	MatrixXd Jacobian = MatrixXd::Zero(6, angles.size() / 2);
	
	for(double t = 0.0; t < TMAX + twidth; t += twidth)
	{
		angles << M_PI / 8 * sin(t *100), M_PI / 8 * sin(t *100),
		          M_PI / 8 * sin(t *100), M_PI / 8 * sin(t *100),
				  M_PI / 8 * sin(t *100), M_PI / 8 * sin(t *100),
		          M_PI / 8 * sin(t *100), M_PI / 8 * sin(t *100);
		
		dangles = angles - angles_past;
		angles_past = angles;
		
		x_Lfoot = TRANS(X.forwardKine(-1.0, angles));
		x_Rfoot = TRANS(X.forwardKine(1.0, angles));
		
		dx_Lfoot = x_Lfoot - x_past_Lfoot;
		dx_Rfoot = x_Rfoot - x_past_Rfoot;
		
		x_past_Lfoot = x_Lfoot;
		x_past_Rfoot = x_Rfoot;
		
		Jacobian = X.getJacobian(-1.0, angles);
		cout << (Jacobian * dangles.block<4,1>(0,0)).transpose() << endl << dx_Lfoot.transpose() << endl;
		
		Jacobian = X.getJacobian(1.0, angles);
		cout << (Jacobian * dangles.block<4,1>(4,0)).transpose() << endl << dx_Rfoot.transpose() << endl;
		
//		cout << X.forwardKine(angles) << endl;
		
		COMJacobian = X.getCOMJacobian(angles);
//		cout << "COMJacobian" << endl << COMJacobian << endl;
		
		COMx		 = X.COMpos(angles);
		
//		dCOMx		 = (COMx - COMx_past) / twidth;
//		COMx_past	 = COMx;
		
//		dangles		 = (angles - angles_past) / twidth;
//		angles_past	 = angles;
		
//		dcomx = COMJacobian * dangles;
		
//		error = dCOMx.head(3) - dcomx.head(3);
		
/*		
		cout << "dangles  " << dangles.transpose() << endl;
		cout << "dCOMx  " << dCOMx.transpose() << endl;
		cout << "dcomx  " << dcomx.transpose() << endl;
		cout << "error  " << error.transpose() << endl;
*/		
//		COMJacobian = X.testGetJacobian(dangles, dCOMx.head(3));	
	}
	
//	for (int i = 0; i < angles.size() + 1; i ++)	cout << "forwardKine(" << i << ")---------------------" << endl <<X.forwardKine(angles,i) << endl;
	return 0;
}
