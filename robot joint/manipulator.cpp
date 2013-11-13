#include <iostream>
#include "manipulator.h"
#include "constant.h"
//#include "Link.h"
#include "Eigen/Dense"

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


//Low Pass Filter
void baseClass::LPF(Vector6d forceDis)
{
	Vector6d dDis = (forceDis-Dis)*W;	
	Dis += dDis*twidth;
}

void baseClass::testJacobian(Vector3d v, VectorXd omega, MatrixXd J)
{
	Vector3d v2, err;
	v2 = (J * omega).head(3);
	err = v - v2;
	cout << v.transpose() << endl << v2.transpose() << endl << err.transpose() << "	<--------- error" << endl << endl;
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
	int dof = 8;
	VectorXd l(dof);
	l << 		L1,		  0.0,	 L2,	 L3,	   -R1,	      0.0,	 R2,	 R3;
	VectorXd m(dof);
	m << -M_PI / 2,  M_PI / 2,	0.0,	0.0, -M_PI / 2,  M_PI / 2,	0.0,	0.0;
	VectorXd n(dof);
	n <<	   0.0,		  0.0,	0.0,	0.0,	   0.0,	      0.0,	0.0,	0.0;
	VectorXd o(dof);
	o << 	   0.0,  M_PI / 2,	0.0,	0.0,       0.0,  M_PI / 2,	0.0,	0.0;
	VectorXd p(dof);
	p << 		M1,		  0.0,    M2,	 M3,  	    M1,	      0.0,	 M2,	 M3;
	
	for(int i = 0; i < dof; i++) links.push_back(new Link(l(i), m(i), n(i), o(i), p(i)));
}
/*
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
*/
Matrix4d manipulator::forwardKine(double LorR, VectorXd angles)
{
	if (LorR < 0)	return forwardKine(angles, links.size() / 2);		//<---Left foot
	else 			return forwardKine(angles, links.size());			//<---Right foot
}


Matrix4d manipulator::forwardKine(VectorXd angles, int to_idx)
{
	Matrix4d T = Matrix4d::Identity();
	int z;
	if (to_idx <= links.size() / 2) z = 0;
	else z = links.size() / 2;
	
	for(int i = z; i < to_idx; i++)	
	{
		assert (i <= links.size());
		T *= links[i]->Transform(angles(i)); 
	}
	
	return T;
}

//get Jacobian
MatrixXd manipulator::getJacobian(double LorR, VectorXd angles)
{
	MatrixXd Jacobian(6, angles.size() / 2);
	Vector3d z, p, p_minus1, pn;
	int first_idx;
	Matrix4d T = Matrix4d::Zero();
	
	if(LorR < 0)
	{
		first_idx = 0;
		pn = TRANS(forwardKine(-1, angles));
	}
	else
	{
		first_idx = angles.size() / 2;
		pn = TRANS(forwardKine(1, angles));
	}
	
	for(int i = 0; i < angles.size() / 2; i++)
	{
		if (i == 0)		T = forwardKine(angles, i);
		else 			T = forwardKine(angles, i + first_idx);
		
		p_minus1 = TRANS(T);
		z = ROT(T).col(2);
		
		p = pn - p_minus1;
		
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
VectorXd manipulator::moment(VectorXd angles, Vector3d fRef, int LorR)
{
	VectorXd n = VectorXd::Zero(3);
	Vector3d x;
	
	for(int i = 0; i < angles.size(); i++)
	{
		if(LorR == 0)	x = TRANS(forwardKine(angles,angles.size() / 2));
		else 			x = TRANS(forwardKine(angles,angles.size()));
		n += x.cross(fRef);
	}
	return n;
}
/*
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
*/

//Center of Mass Position calculation
Vector4d manipulator::COMpos(int start, VectorXd angles)		//start -> joint number
{
	Vector4d ans = Vector4d::Zero();
	Vector4d c[angles.size()]; 
	c[0] << -L1/2, 0.0, 0.0, 1.0;
	c[1] <<   0.0, 0.0, 0.0, 1.0;
	c[2] << -L2/2, 0.0, 0.0, 1.0;
	c[3] << -L3/2, 0.0, 0.0, 1.0;
	c[4] <<  R1/2, 0.0, 0.0, 1.0;
	c[5] <<   0.0, 0.0, 0.0, 1.0;
	c[6] << -R2/2, 0.0, 0.0, 1.0;
	c[7] << -R3/2, 0.0, 0.0, 1.0;
	
	double M_total = 0.0;
	int terminal;
	if (start < links.size() / 2)	terminal = links.size() / 2; 
	else 							terminal = links.size();
	
    for (int h = start; h < terminal; h ++)
    {
        M_total += links[h]->mass;
    }
    
	for(int i = start; i < terminal; i++)
	{
		c[i] = forwardKine(angles, i+1) * c[i];
//		c[i].head(3) += forwardKine(angles, i).block<3,1>(0,3);
	}
	
	for(int j = start; j < terminal; j++)
	{
		ans += (c[j] * links[j]->mass) / M_total;
	}
	ans(3) = M_total;
	
	return ans;
}

Vector4d manipulator::COMpos(VectorXd angles)
{
	Vector4d ans = Vector4d::Zero();
	Vector4d Lcom = COMpos(0, angles);			//<---- Lfoot COM
	Vector4d Rcom = COMpos(4, angles);			//<---- Rfoot COM
	double M_total = Lcom(3) + Rcom(3);
	ans.head(3) = (Lcom.head(3) * Lcom(3) + Rcom.head(3) * Rcom(3)) / M_total;
	ans(3) = M_total;
	
	return ans;
}

//COM Jacobian
MatrixXd manipulator::getCOMJacobian(VectorXd angles)
{
	MatrixXd COMJacobian = MatrixXd::Zero(6, angles.size());
	Matrix4d T = Matrix4d::Zero();
	Vector4d COM[angles.size()];
	Vector3d z, p_minus, p;
	double M_total = 0.0;
	
	for(int i = 0; i < angles.size(); i++)
	{
		COM[i] = COMpos(i, angles);
		M_total += links[i]->mass;
	}
	
	for(int i = 0; i < angles.size(); i++)
	{
		if(i < angles.size() / 2)			T = forwardKine(angles, i);
		else if(i == angles.size() / 2) 	T = forwardKine(angles, 0);
		else 								T = forwardKine(angles, i);
		p_minus = TRANS(T);
		p = COM[i].head(3) - p_minus;
		
		z = ROT(T).col(2);
		
		COMJacobian.block<3,1>(0,i) = (COM[i](3) / M_total) * z.cross(p);
		COMJacobian.block<3,1>(3,i) = z; 
	}
	
	return COMJacobian;
}

MatrixXd manipulator::getLfoot2COMJacobian(VectorXd angles)
{
	Matrix4d T = forwardKine(-1.0, angles);
	Matrix3d L2Brot = ROT(T).transpose();
	Vector3d LfootPos = TRANS(T);
	Vector4d COMPos = COMpos(angles);
	
	Vector3d L2C_base = COMPos.head(3) - LfootPos;
	Vector3d L2C_Lfoot = L2Brot * L2C_base;
	Vector3d z;
	
	MatrixXd JacL2C = MatrixXd::Zero(6, angles.size());
	MatrixXd JCOM = getCOMJacobian(angles);
	MatrixXd JKine = MatrixXd::Zero(6, angles.size());
	
	MatrixXd JKine_L = getJacobian(-1.0, angles);
	MatrixXd JKine_R = getJacobian(1.0, angles);
	
    for(int i = 0; i < angles.size(); i++)
	{
		if(i < angles.size() / 2)	JKine.block<3, 1>(0, i) = JKine_L.block<3, 1>(0, i);
		else 						JKine.block<3, 1>(0, i) = JKine_L.block<3, 1>(0, i - angles.size() / 2);
		z = L2Brot * JKine.block<3, 1>(3, i);
		if(i < angles.size ()/ 2)	JacL2C.block<3, 1>(0, i) = L2Brot * (JCOM.block<3, 1>(0, i) - JKine.block<3, 1>(0, i)) - z.cross(L2C_Lfoot);
		else 						JacL2C.block<3, 1>(0, i) = L2Brot * (JCOM.block<3, 1>(0, i)) - z.cross(L2C_Lfoot);
		JacL2C.block<3, 1>(3, i) = z;
    }
//	cout << "JCOM" << endl << JCOM << endl << "JKine" << endl << JKine << endl;
	return JacL2C;
}

MatrixXd manipulator::getLfoot2RfootJacobian(VectorXd angles)
{
	MatrixXd Jac_RinB = MatrixXd::Zero(6, angles.size());
	MatrixXd Jac_RinL = MatrixXd::Zero(6, angles.size());
	MatrixXd Convert = MatrixXd::Zero(6, 6);
	Matrix3d Rot_B2L = ROT(forwardKine(-1.0, angles)).transpose();
	
	Convert.block<3,3>(0,0) = Rot_B2L;
	Convert.block<3,3>(3,3) = Rot_B2L;
	
	Jac_RinB.block<6,4>(0,0) = -getJacobian(-1.0, angles);
	Jac_RinB.block<6,4>(0,4) = getJacobian(1.0, angles);
	Jac_RinL = Convert * Jac_RinB;
	
	return Jac_RinL;
}
