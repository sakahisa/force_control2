#include <iostream>
#include <cmath>
#include "Eigen/Dense"
#include <vector>
#include "Set.h"

#define twidth 0.00001	//sample time
#define TMAX 0.001
#define l1 1.0			//length of Link1
#define l2 0.5			//length of Link2
#define l3 0.5			//length of Link3
#define W 1000			//cut off
#define M 1.0			//mass
#define D 2.5			//coefficient of viscosity

#define TRANS(x) x.block<3,1>(0,3)
#define ROT(x) x.block<3,3>(0,0) 

using namespace std;
using namespace Eigen;

template<typename _Matrix_Type_>
//SVD
bool svdInverse(const _Matrix_Type_ &a, _Matrix_Type_ &result, double epsilon = std::numeric_limits<typename _Matrix_Type_::Scalar>::epsilon())
{
	if(a.rows()<a.cols())
		return false;
	
	Eigen::JacobiSVD< _Matrix_Type_ > svd = a.jacobiSvd(Eigen::ComputeFullU |Eigen::ComputeFullV);
	
	typename _Matrix_Type_::Scalar tolerance = epsilon * std::max(a.cols(),a.rows()) * svd.singularValues().array().abs().maxCoeff();
	
	result = svd.matrixV() * _Matrix_Type_( (svd.singularValues().array().abs() >
	tolerance).select(svd.singularValues().array().inverse(), 0) ).asDiagonal() * svd.matrixU().adjoint();
}

typedef Eigen::Matrix<double, 6, 1> Vector6d;

//Low Pass Filter
void baseClass::LPF(Vector6d forceDis)
{
    Vector6d dDis = (forceDis-Dis)*W;
    Dis += dDis*twidth;
}

Link::Link(double l,double m,double n,double o, double p)
{
	a = l;
	alpha = m;
	offset = n;
	theta_offset = o;
    mass = p;
}
//Transform Matrix
Matrix4d Link::Transform(double angle)
{
	Matrix4d T;
	
	double theta = angle + theta_offset;
	
/*
	T(0,0)=cos(theta);					T(0,1)=-sin(theta);				T(0,2)=0.0;				T(0,3)=a;
	T(1,0)=sin(theta)*cos(alpha);		T(1,1)=cos(theta)*cos(alpha);	T(1,2)=-sin(alpha);		T(1,3)=-sin(alpha)*offset;
	T(2,0)=sin(theta)*sin(alpha);		T(2,1)=cos(theta)*sin(alpha);	T(2,2)=cos(alpha);		T(2,3)=cos(alpha)*offset;
	T(3,0)=0.0;							T(3,1)=0.0;						T(3,2)=0.0;				T(3,3)=1.0;
*/	
	
	T(0,0) = cos(theta);	T(0,1) = -sin(theta)*cos(alpha);		T(0,2) = sin(theta)*sin(alpha);		T(0,3) = a*cos(theta);
	T(1,0) = sin(theta);	T(1,1) = cos(theta)*cos(alpha);		T(1,2) = -cos(theta)*sin(alpha);		T(1,3) = a*sin(theta);
	T(2,0) = 0.0;		T(2,1) = sin(alpha);			T(2,2) = cos(alpha);			T(2,3) = offset;
	T(3,0) = 0.0;		T(3,1) = 0.0;				T(3,2) = 0.0;				T(3,3) = 1.0;
	return T;
}

force::force()
{
	Vector6d Z = VectorXd::Zero(6);
	Cmd = Z;							//command
	Res = Z;							//response
	Dis = Z;							//disturbance
	dDis = Z;							//difference of disturbance
}

tau::tau()
{
	VectorXd Z = VectorXd::Zero(3);
	Cmd = Z;
	Res = Z;
	Dis = Z;
}

//viscous friction
void tau::Disturbance(VectorXd omega){
	Dis=-D*omega;
}

manipulator::manipulator()
{
	int dof = 3;
	VectorXd l(3);
	l << l1, l2, l3;
	VectorXd m(3);
	m << 0.0, 0.0, 0.0;
	VectorXd n(3);
	n << 0.0, 0.0, 0.0;
	VectorXd o(3);
	o << 0.0, 0.0, 0.0;
    VectorXd p(3);
    p << M, M, M;
	
    for(int i = 0; i < dof; i++) links.push_back(new Link(l(i), m(i), n(i), o(i), p(i)));
}

//torque to force
Vector6d manipulator::torque2Force(VectorXd angles, VectorXd tauRef)
{
	return pseudoInverse(getJacobian(angles)).transpose()*tauRef;
}

//force to torque
VectorXd manipulator::force2Torque(VectorXd angles, Vector6d forceRef)
{
	return getJacobian(angles).transpose()*forceRef;
}

//force controll
VectorXd manipulator::forceControl(VectorXd angles, Vector6d forceRef, VectorXd omega)
{
	Vector6d forceDis;
	
	T.Disturbance(omega);
	
	F.Cmd = forceRef + F.Dis;
	T.Cmd = force2Torque(angles,F.Cmd);
	T.Res = T.Cmd+T.Dis;
	F.Res = torque2Force(angles,T.Res);		//sensor force
	
	forceDis = F.Cmd-F.Res;
	F.LPF(forceDis);
	
	return F.Res;
}

Matrix4d manipulator::forwardKine(VectorXd angles)
{
	return forwardKine(angles, links.size());
}

Matrix4d manipulator::forwardKine(VectorXd angles, int to_idx)
{
	Matrix4d T = Matrix4d::Identity();
	
	for(int i = 0; i < to_idx; i++)	
	{
		assert (i <= links.size());
		T *= links[i]->Transform(angles(i)); 
	}

	return T;
}

//get Jacobian
MatrixXd manipulator::getJacobian(VectorXd angles)
{
	MatrixXd Jacobian(6, angles.size());
	Vector3d z, p, p_minus1, pn;
	
	pn = TRANS(forwardKine(angles));
	
	for(int i = 0; i < links.size(); i++)
	{
		p_minus1 = TRANS(forwardKine(angles, i));
		p = pn - p_minus1;
		
		z = ROT(forwardKine(angles, i)).col(2);
		
		Jacobian.block<3,1>(0,i) = z.cross(p);
		Jacobian.block<3,1>(3,i) = z; 
	}
	return Jacobian;
}

MatrixXd manipulator::pseudoInverse(MatrixXd Jacobian)
{
	MatrixXd Inv;
	MatrixXd JtJ = Jacobian.transpose()*Jacobian;
	svdInverse(JtJ, Inv);
	return Inv*Jacobian.transpose();	
}

//calculation moment
VectorXd manipulator::moment(VectorXd angles, Vector3d fRef)
{
	VectorXd n = VectorXd::Zero(3);
	Vector3d x;
	
	for(int i = 0; i < angles.size(); i++)
	{
		x = TRANS(forwardKine(angles,3));
		n += x.cross(fRef);
	}
	
	return n;
}

VectorXd manipulator::invKine(Vector3d Error, VectorXd angles)
{
	double K = Error.norm()*1000;
//	cout << K << endl;
	if(Error.norm() < 0.0001)			//stop
	{
		return angles;
	}
	
	
	angles += twidth * K * pseudoInverse(ROT(getJacobian(angles))) * Error;
	return angles;
}

//Center of Mass Position calculating
Vector4d manipulator::COMpos(int start, VectorXd angles)		//start -> joint number
{
	Link A;
	VectorXd ans = VectorXd::Zero(4);
	Vector4d c[3]; 
	c[0] << -l1/2, 0.0, 0.0, 1.0;
	c[1] << -l2/2, 0.0, 0.0, 1.0;
	c[2] << -l3/2, 0.0, 0.0, 1.0;
	
    double M_total = 0.0;
    for (int h = start; h < links.size(); h ++)
    {
        M_total += links[h]->mass;
    }
	
	for(int i = start; i < 3; i++)
	{
		c[i] = forwardKine(angles, i+1) * c[i];
//		c[i].head(3) += forwardKine(angles, i).block<3,1>(0,3);
	}
	for(int j = start; j < 3; j++)
	{
		ans += (c[j] * M) / M_total;
	}
	ans(3) = M_total;
	
//	cout << c[0] << endl;
	return ans;
}

//COM Jacobian
MatrixXd manipulator::getCOMJacobian(VectorXd angles)
{
    MatrixXd COMJacobian(6, angles.size());
    Vector4d COM[angles.size()];
    Vector3d z, p_minus, p;

    for(int i = 0; i < angles.size(); i++)
    {
        COM[i] = COMpos(i, angles);
    }

    for(int i = 0; i < links.size(); i++)
    {
        p_minus = TRANS(forwardKine(angles, i));
        p = COM[i].head(3) - p_minus;

        z = ROT(forwardKine(angles, i)).col(2);

        COMJacobian.block<3,1>(0,i) = (COM[i](3) / (M * 2)) * z.cross(p);
        COMJacobian.block<3,1>(3,i) = z;
    }

    return COMJacobian;
}
