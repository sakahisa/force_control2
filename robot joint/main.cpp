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
	
	Vector3d v_L2R;
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
	MatrixXd Jacobian_L = MatrixXd::Zero(6, angles.size() / 2);
	MatrixXd Jacobian_R = MatrixXd::Zero(6, angles.size() / 2);
	for(double t = 0.0; t < TMAX + twidth; t += twidth)
	{
		angles << M_PI / 8 * sin(t), M_PI / 8 * sin(t),
		          M_PI / 8 * sin(t), M_PI / 8 * sin(t),
				  M_PI / 8 * sin(t), M_PI / 8 * sin(t),
		          M_PI / 8 * sin(t), M_PI / 8 * sin(t);
		
		dangles = angles - angles_past;
		angles_past = angles;
		
		x_Lfoot = TRANS(X.forwardKine(-1.0, angles));
		x_Rfoot = TRANS(X.forwardKine(1.0, angles));
		COMx = (X.COMpos(angles)).head(3);
		
		dx_Lfoot = x_Lfoot - x_past_Lfoot;
		dx_Rfoot = x_Rfoot - x_past_Rfoot;
		dCOMx = COMx - COMx_past;
		
		x_past_Lfoot = x_Lfoot;
		x_past_Rfoot = x_Rfoot;
		COMx_past = COMx;
		
		Jacobian_L = X.getJacobian(-1.0, angles);
		
		Jacobian_R = X.getJacobian(1.0, angles);
		
		COMJacobian = X.getCOMJacobian(angles);
		
		v_L2R = (ROT(X.forwardKine(-1.0, angles))).transpose() * (dx_Rfoot - dx_Lfoot);
		v_COM = (ROT(X.forwardKine(-1.0, angles))).transpose() * (dCOMx / twidth);
		
		cout << "v_L2R					v_COM" << endl;
		cout << v_L2R.transpose() << "		" << v_COM.transpose() << endl << endl;
		
		v_L2R = (X.getLfoot2RfootJacobian(angles) * dangles).head(3);
		v_COM = (X.getLfoot2COMJacobian(angles) * dangles).head(3);
		
		cout << v_L2R.transpose() << "		" << v_COM.transpose() << endl << endl;
		
		
//		cout << X.forwardKine(angles) << endl;
		
//		cout << "COMJacobian" << endl << COMJacobian << endl;
		
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
