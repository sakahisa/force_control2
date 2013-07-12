#include <iostream>
#include <cmath>
#include <fstream>
#include "Eigen/Dense"
#include <vector>

#define twidth 0.001
#define TMAX 5.0
#define l1 1.0
#define l2 0.5
#define l3 0.5
#define W 500		//cut off
#define M 1.0


#define TRANS(x) x.block<3,1>(0,3)
#define ROT(x) x.block<3,3>(0,0) 

using namespace std;
using namespace Eigen;

typedef Eigen::Matrix<double, 6, 1> Vector6d;

class Link
{
	public:
	double a,alpha,offset,theta_offset;
	
	Link(){};
	Link(double l,double m,double n,double o);
	Matrix4d Transform(double angle);
};

class force
{
	public:
	Vector6d Cmd,Res,Dis,dDis;
	
	force();
	void LPF(Vector6d forceDis);
};

class tau
{
	public:
	VectorXd Cmd,Res,Dis;
	
	tau();
};

class manipulator
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
	VectorXd forceControl(VectorXd angles, Vector6d forceRef);
	
	MatrixXd pseudoInverse(MatrixXd Jacobian);
	VectorXd moment(VectorXd angles, Vector3d fRef);
	
	std::vector<Link* > links;
};



Link::Link(double l,double m,double n,double o)
{
	a = l;
	alpha = m;
	offset = n;
	theta_offset = o;
}

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

force::force(){
	Vector6d Z = VectorXd::Zero(6);
	Cmd = Z;
	Res = Z;
	Dis = Z;
	dDis = Z;
}

void force::LPF(Vector6d forceDis)
{
	Dis = forceDis;
	
	dDis = (forceDis-Dis)*W;
	
	Dis += dDis*twidth;
}

tau::tau(){
	VectorXd Z = VectorXd::Zero(3);
	Cmd = Z;
	Res = Z;
	Dis = Z;
}

manipulator::manipulator()
{
	int dof = 3;
	VectorXd l(3);
	l << l1, l2, l3;
	VectorXd m(3);
	m << 0.0, M_PI/2, 0.0;
	VectorXd n(3);
	n << 0.0, 0.0, 0.0;
	VectorXd o(3);
	o << 0.0, 0.0, 0.0;
	
	for(int i = 0; i < dof; i++) links.push_back(new Link(l(i),m(i),n(i),o(i)));
	
}

Vector6d manipulator::torque2Force(VectorXd angles, VectorXd tauRef)
{
	return pseudoInverse(getJacobian(angles)).transpose()*tauRef;
}
VectorXd manipulator::force2Torque(VectorXd angles, Vector6d forceRef)
{
	return getJacobian(angles).transpose()*forceRef;
}

VectorXd manipulator::forceControl(VectorXd angles, Vector6d forceRef)
{
	/*
	 * Naoki Oda ?
	 */
	Vector6d forceDis;
	
	F.Cmd = forceRef + F.Dis;
	T.Cmd = force2Torque(angles,F.Cmd);
	T.Res = T.Cmd+T.Dis;
	F.Res = torque2Force(angles,T.Res);		//sensor force
	
	forceDis = F.Cmd-F.Res;
	F.LPF(forceDis);
	
	return F.Res;
}



Matrix4d manipulator::forwardKine(VectorXd angle)
{
	return forwardKine(angle, links.size());
}

Matrix4d manipulator::forwardKine(VectorXd angle, int to_idx)
{
	Matrix4d T = Matrix4d::Identity();
	
	for(int i = 0; i < to_idx; i++)	
	{
		assert (i <= links.size());
		T *= links[i]->Transform(angle(i) ); 
	}

	return T;
}

MatrixXd manipulator::getJacobian(VectorXd angle)
{
	MatrixXd Jacobian(6, angle.size());
	Vector3d z, p, p_minus1, pn;
	
	pn = TRANS(forwardKine(angle));
	
	for(int i = 0; i < links.size(); i++)
	{
		p_minus1 = TRANS(forwardKine(angle, i));
		p = pn - p_minus1;
		
		z = ROT(forwardKine(angle, i)).col(2);
		
		Jacobian.block<3,1>(0,i) = z.cross(p);
		Jacobian.block<3,1>(3,i) = z; 
	}
	return Jacobian;
}

MatrixXd manipulator::testGetJacobian(VectorXd dangles, Vector3d dx)
{
	MatrixXd Jacobian(6,dangles.size());
	
	for(int i = 0; i<3; i++)
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
	return (Jacobian.transpose()*Jacobian).inverse()*Jacobian.transpose();
}

VectorXd manipulator::moment(VectorXd angles, Vector3d fRef)
{
	VectorXd n = VectorXd::Zero(3);
	Vector3d x;
	
	for(int i=0; i<angles.size(); i++)
	{
		x = TRANS(forwardKine(angles,3));
		n += x.cross(fRef);
	}
	
	return n;
} 

int main()
{
	ofstream ofs_f("force.txt");
	ofstream ofs_p("position.txt");
	
	double t;
	Link A;
	manipulator X;	
	
	VectorXd angles(3);
//	VectorXd angles2 = VectorXd::Zero(3);
//	VectorXd dangles(3);
	
	angles << M_PI/6, M_PI/4, M_PI/3;
	
//	MatrixXd J(6,angles.size());
//	J=X.getJacobian(angles);
	
	VectorXd omega(3);

//	MatrixXd J(6,angles.size());
//	MatrixXd Jin(angles.size(),6);
	
	Vector3d x = TRANS(X.forwardKine(angles,3)); 
//	VectorXd x2 = TRANS(X.forwardKine(angles,3));
//	Vector3d dx;// std::vector< Vector3d >
	VectorXd a = VectorXd::Zero(3);
	VectorXd v = VectorXd::Zero(3);	
	
	Vector6d forceRef,forceRes;
	
	Vector3d fRef;
	fRef << 1.0, 1.0, 1.0;
	
	forceRef.block<3,1>(0,0) = fRef;
	
//	forceRef.block<3,1>(3,0) = X.moment(angles,fRef);
	
//	forceRes = X.forceControl(angles,forceRef);
	
//	cout << forceRef.transpose() << endl << forceRes.transpose() << endl << endl;
	

	for(t = 0.0; t < TMAX; t += twidth){
		x = TRANS(X.forwardKine(angles,3)); 
		forceRef.block<3,1>(3,0) = X.moment(angles,fRef);
		
		forceRes = X.forceControl(angles,forceRef);
		
		ofs_f << t << " " << forceRef.transpose() << " " << forceRes.transpose() << endl;
		ofs_p << t << " " << x.transpose() << endl;
		
		a = forceRes.head(3)/M;
		v += a*twidth;
		omega = X.pseudoInverse(ROT(X.getJacobian(angles)))*v;
		
		angles += omega*twidth;

/*
		dx = (x-x2)/twidth;
		x2 = x;
		
		dangles = (angles-angles2)/twidth;
		angles2 = angles;
		
		J = X.getJacobian(angles);
		v = J*dangles;
		Jin = X.pseudoInverse(J);
		omega = ROT(J).inverse()*v.head(3);

		cout << "confirm Jacobian Effect ----> v:dx" <<endl<< v.head(3).transpose() <<endl<<  dx.transpose()  << endl;		//confirm Jacobian Effect
		cout << "confirm Jacobian.inverse Effect ----> dangles,omega" <<endl<< dangles.transpose() <<endl<<  omega.transpose()  << endl;		//confirm Jacobian Effect
*/
	}
	
	return 0;
}
