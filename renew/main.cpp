#include <iostream>
#include <cmath>
#include <fstream>
#include "Eigen/Dense"
#include <vector>

#define twidth 	0.000001	//sample time
#define TMAX 	0.01
#define l1 		1.0			//length of Link1
#define l2 		1.0			//length of Link2
#define l3 		1.0			//length of Link3
#define W 		1000			//cut off
#define M 		1.0			//mass
#define D 		2.5			//coefficient of viscosity

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
VectorXd Zero5 = VectorXd::Zero(5);

class baseClass
{
	public:
	Vector6d Dis;
	void LPF(Vector6d forceDis);
};

//Low Pass Filter
void baseClass::LPF(Vector6d forceDis)
{
	Vector6d dDis = (forceDis-Dis)*W;	
	Dis += dDis*twidth;
}

class Link
{
	public:
	double a,alpha,offset,theta_offset;
	
	Link(){};
	Link(double l,double m,double n,double o);
	Matrix4d Transform(double angle);
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
		
	Matrix4d forwardKine(VectorXd angle);
	Matrix4d forwardKine(VectorXd angle, int to_idx);
	
	MatrixXd getJacobian(VectorXd angle);
	MatrixXd testGetJacobian(VectorXd angle, Vector3d v);
	
	Vector6d torque2Force(VectorXd angles, VectorXd tauRef);
	VectorXd force2Torque(VectorXd angles, Vector6d forceRef);
	VectorXd forceControl(VectorXd angles, Vector6d forceRef, VectorXd omega);
	
	MatrixXd pseudoInverse(MatrixXd Jacobian);
	VectorXd moment(VectorXd angles, Vector3d fRef);
	
	VectorXd invKine(Vector3d Error, VectorXd angles);
	
	VectorXd COMpos(int start, VectorXd angles);
	VectorXd COMpos(VectorXd angles);
	MatrixXd getCOMJacobian(VectorXd angles);
	
	std::vector<Link* > links;
};



Link::Link(double l,double m,double n,double o)
{
	a = l;
	alpha = m;
	offset = n;
	theta_offset = o;
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
	
	for(int i = 0; i < dof; i++) links.push_back(new Link(l(i),m(i),n(i),o(i)));
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

MatrixXd manipulator::testGetJacobian(VectorXd dangles, Vector3d dx)
{
	MatrixXd Jacobian(6,dangles.size());
	
	for(int i = 0; i<dx.size(); i++)
	{
		for(int j = 0; j<dangles.size(); j++)
		{
			Jacobian(i,j) = dx(i) / dangles(j);
		}
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

//Center of Mass Position calculation
VectorXd manipulator::COMpos(int start, VectorXd angles)		//start -> joint number
{
	VectorXd ans = Zero5;
	Vector4d c[angles.size()]; 
	c[0] << -l1/2, 0.0, 0.0, 1.0;
	c[1] << -l2/2, 0.0, 0.0, 1.0;
	c[2] << -l3/2, 0.0, 0.0, 1.0;
	
	double M_total = M * (angles.size() - start);
	
	for(int i = start; i < angles.size(); i++)
	{
		c[i] = forwardKine(angles, i+1) * c[i];
//		c[i].head(3) += forwardKine(angles, i).block<3,1>(0,3);
	}
	for(int j = start; j < angles.size(); j++)
	{
		ans.head(3) += (c[j].head(3) * M) / M_total;
	}
	ans(3) = 1;
	ans(4) = M_total;
	
	return ans;
}

VectorXd manipulator::COMpos(VectorXd angles)
{
	return COMpos(0, angles);
}

//COM Jacobian
MatrixXd manipulator::getCOMJacobian(VectorXd angles)
{
	MatrixXd COMJacobian(6, angles.size());
	VectorXd COM[angles.size()];
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
/*		
		cout << "i       " << i << endl;
		cout << "COM[" << i << "]  " << COM[i].transpose() << endl;
		cout << "p       " << p.transpose() << endl;
		cout << "p_minus " << p_minus.transpose() << endl;
		cout << "z       " << z.transpose() << endl;
*/		
		COMJacobian.block<3,1>(0,i) = (COM[i](4) / (M * 3)) * z.cross(p);
		COMJacobian.block<3,1>(3,i) = z; 
	}
	
	return COMJacobian;
}

int main()
{
	Link A;
	manipulator X;	
	
	ofstream ofs("edgepos.txt");
	
	Vector3d angles = Vector3d::Zero();
	
	Vector3d x = TRANS(X.forwardKine(angles,3)); 
	VectorXd v = VectorXd::Zero(3);
	VectorXd a = VectorXd::Zero(3);
	
	Vector3d xRef, Error;
	
	VectorXd dangles;
	VectorXd angles_past = VectorXd::Zero(3);
	VectorXd dCOMx = Zero5;
	VectorXd COMx = Zero5;
	VectorXd COMx_past = Zero5;
	
	
	MatrixXd COMJacobian = MatrixXd::Zero(6, angles.size());
	
	for(double t = 0.0; t < TMAX; t += twidth)
	{
		angles << M_PI / 4 * sin(t *1000), -M_PI / 4 * sin(t *1000), -M_PI / 4 * sin(t *1000);
/*		
		Error = xRef - x;
		angles = X.invKine(Error, angles);
		x = TRANS(X.forwardKine(angles,3));
*/		
//		cout << X.forwardKine(angles) << endl;
		
		COMJacobian = X.getCOMJacobian(angles);
//		cout << COMJacobian << endl;
		
		COMx = X.COMpos(angles);
		
		dangles = angles - angles_past;
		dCOMx = COMx - COMx_past;
		COMJacobian = X.testGetJacobian(dangles, dCOMx.head(3));
		angles_past = angles;
		COMx_past = COMx;
//		cout << "----------------test-------------------" << endl << COMJacobian << endl << "---------------------------------------" << endl;
		ofs << t << " " << COMx.transpose() << " " << TRANS(X.forwardKine(angles)).transpose() << endl;
	}
/*	
	for(int i = 0; i < 5; i ++)
	{
		cout << X.forwardKine(angles, i) << endl << endl;
	}
	*/
	return 0;
}
