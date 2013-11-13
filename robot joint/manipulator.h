#include "Eigen/Dense"
#include "Link.h"
#include <vector>

using namespace Eigen;

typedef Eigen::Matrix<double, 6, 1> Vector6d;

class baseClass
{
	public:
	Vector6d Dis;
	void LPF(Vector6d forceDis);
	void testJacobian(Vector3d v, VectorXd omega, MatrixXd J);
};

class force : public baseClass
{
	public:
	Vector6d Cmd,Res,Dis,dDis;
	
	force();
};

class tau
{
	public:
	VectorXd Cmd,Res,Dis;
	
	tau();
	void Disturbance(VectorXd omega);
};


class manipulator : public baseClass
{
	public:
	
	force F;
	tau T;
	
	manipulator();
		
	Matrix4d forwardKine(double LorR, VectorXd angle);
	Matrix4d forwardKine(VectorXd angle, int to_idx);
	
	MatrixXd getJacobian(double LorR, VectorXd angle);
	
	Vector6d torque2Force(VectorXd angles, VectorXd tauRef);
	VectorXd force2Torque(VectorXd angles, Vector6d forceRef);
	VectorXd forceControl(VectorXd angles, Vector6d forceRef, VectorXd omega);
	
	MatrixXd pseudoInverse(MatrixXd Jacobian);
	VectorXd moment(VectorXd angles, Vector3d fRef, int LorR);
	
	VectorXd invKine(Vector3d Error, VectorXd angles);
	
	Vector4d COMpos(int start, VectorXd angles);
	Vector4d COMpos(VectorXd angles);
	MatrixXd getCOMJacobian(VectorXd angles);
	
	MatrixXd getLfoot2COMJacobian(VectorXd angles);
	MatrixXd getLfoot2RfootJacobian(VectorXd angles);
	
	std::vector<Link* > links;
};
