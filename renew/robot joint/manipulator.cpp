#include <iostream>
#include "manipulator.h"
#include "constant.h"
//#include "Link.h"
#include "Eigen/Dense"

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
		if (i == 0)
		{
			p_minus1 = TRANS(forwardKine(angles, i));
			z = ROT(forwardKine(angles, i)).col(2);
		}
		else
		{
			p_minus1 = TRANS(forwardKine(angles, i + first_idx));
			z = ROT(forwardKine(angles, i + first_idx)).col(2);
		}
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
	Vector4d Lcom = COMpos(0, angles);
	Vector4d Rcom = COMpos(4, angles);
	double M_total = Lcom(3) + Rcom(3);
	ans.head(3) = (Lcom.head(3) * Lcom(3) + Rcom.head(3) * Rcom(3)) / M_total;
	ans(3) = M_total;
	
	std::cout << "Lcom	" << Lcom.transpose() << std::endl;
	std::cout << "Rcom	" << Rcom.transpose() << std::endl;
	std::cout << "COMp	" << ans.transpose()  << std::endl;
	std::cout << "links.size() " << links.size() << std::endl;
	
	return ans;
}

//COM Jacobian
MatrixXd manipulator::getCOMJacobian(VectorXd angles)
{
	MatrixXd COMJacobian = MatrixXd::Zero(6, angles.size());
	Vector4d COM[angles.size()];
	Vector3d z, p_minus, p;
	double M_total;
	
	for(int i = 0; i < angles.size(); i++)
	{
		COM[i] = COMpos(i, angles);
		M_total += links[i]->mass;
	}
	
	for(int i = 0; i < links.size(); i++)
	{
		p_minus = TRANS(forwardKine(angles, i));
		p = COM[i].head(3) - p_minus;
		
		z = ROT(forwardKine(angles, i)).col(2);
/*		
		cout << "i       " << i << endl;
		cout << "COM[" << i << "]  " << COM[i].transpose() << endl;
		cout << "p       " << p.transpose() << endl;
		cout << "p_minus " << p_minus.transpose() << endl;
		cout << "z       " << z.transpose() << endl;
*/		
		if(i < links.size()-1) COMJacobian.block<3,1>(0,i) = (COM[i](3) / M_total) * z.cross(p);
		COMJacobian.block<3,1>(3,i) = z; 
	}
	
	return COMJacobian;
}

MatrixXd manipulator::getLfoot2COMJacobian(VectorXd angles)
{
	Vector3d LfootPos = TRANS(forwardKine(-1.0, angles));
	Vector4d COMPos = COMpos(angles);
	
	
}


MatrixXd manipulator::getLfoot2RfootJacobian(VectorXd angles)
{
	MatrixXd Jac_B2R = MatrixXd::Zero(6, angles.size());
	MatrixXd Jac_L2R = MatrixXd::Zero(6, angles.size());
	MatrixXd Exchange = MatrixXd::Zero(6, 6);
	Matrix3d Rot_B2L = ROT(forwardKine(-1.0, angles));
	
	Exchange.block<3,3>(0,0) = Rot_B2L;
	Exchange.block<3,3>(3,3) = Rot_B2L;
	
//	JacB2R = getJacobian(angles, angles.size();
}

Vector3d manipulator::getLfoot2RfootPos(VectorXd angles)
{
	Vector3d LfootPos, RfootPos, L2R_base, L2R_Lfoot;
	Matrix3d Rot_L2R = Matrix3d::Zero();
	
	LfootPos = TRANS(forwardKine(-1.0, angles));
	RfootPos = TRANS(forwardKine(1.0, angles));
	
	L2R_base = LfootPos -RfootPos;
	
	Rot_L2R = (ROT(forwardKine(-1.0, angles))).transpose() * ROT(forwardKine(1.0, angles));
	
	L2R_Lfoot = Rot_L2R * L2R_base;
	
	return L2R_Lfoot;
}
