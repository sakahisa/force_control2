#include <iostream>
#include <cmath>
#include <fstream>
#include "Eigen/Dense"
#include <vector>

#define twidth 0.0001
#define TMAX 0.10
#define l1 1.0
#define l2 0.5
#define l3 0.5

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

class manipulator
{
	public:
	
	manipulator();

	manipulator(VectorXd angle);
		
	Matrix4d forwardKine(VectorXd angle);
	Matrix4d forwardKine(VectorXd angle, int to_idx);
	
	MatrixXd getJacobian(VectorXd angle);
	MatrixXd testGetJacobian(VectorXd angle, Vector3d v);
	
	Vector6d torque2Force(VectorXd tauRef);
	VectorXd force2Torque(Vector6d forceRef);
	VectorXd forceControl(Vector6d forceRef);

	std::vector<Link* > links;
	
};


Link::Link(double l,double m,double n,double o)
{
	a=l;
	alpha=m;
	offset=n;
	theta_offset=o;
}

Matrix4d Link::Transform(double angle)
{
	Matrix4d T;
	
	double theta=angle + theta_offset;
	/*
	T(0,0)=cos(theta);					T(0,1)=-sin(theta);				T(0,2)=0.0;				T(0,3)=a;
	T(1,0)=sin(theta)*cos(alpha);		T(1,1)=cos(theta)*cos(alpha);	T(1,2)=-sin(alpha);		T(1,3)=-sin(alpha)*offset;
	T(2,0)=sin(theta)*sin(alpha);		T(2,1)=cos(theta)*sin(alpha);	T(2,2)=cos(alpha);		T(2,3)=cos(alpha)*offset;
	T(3,0)=0.0;							T(3,1)=0.0;						T(3,2)=0.0;				T(3,3)=1.0;
*/
	T(0,0)=cos(theta);	T(0,1)=-sin(theta)*cos(alpha);		T(0,2)=sin(theta)*sin(alpha);		T(0,3)=a*cos(theta);
	T(1,0)=sin(theta);	T(1,1)=cos(theta)*cos(alpha);		T(1,2)=-cos(theta)*sin(alpha);		T(1,3)=a*sin(theta);
	T(2,0)=0.0;		T(2,1)=sin(alpha);			T(2,2)=cos(alpha);			T(2,3)=offset;
	T(3,0)=0.0;		T(3,1)=0.0;				T(3,2)=0.0;				T(3,3)=1.0;
	return T;
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

manipulator::manipulator(VectorXd angle)
{
	int dof = 3;
	VectorXd l(3);
	l << l1, l2, l3;
	VectorXd m(3);
	m << 0.0, 0.0, 0.0;
	VectorXd n(3);
	n << 0.0, 0.0, 0.0;
	VectorXd o(3);
	o = angle;
	
	for(int i = 0; i < dof; i++) links.push_back(new Link(l(i),m(i),n(i),o(i)));
}

Vector6d manipulator::torque2Force(VectorXd angles, VectorXd tauRef)
{
	return getJacobian(angles).inverse().transpose()*forceRef;
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
	;
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
 
int main()
{
	std::ofstream ofs("position.txt");
	
	int i;
	double t;
	Link A;
	
	VectorXd angles = VectorXd::Zero(3);
	VectorXd angles2(3);
	VectorXd dangles(3);
	
	VectorXd omega;
	VectorXd v(6);
	MatrixXd J(6,angles.size());
	Matrix4d T;	
	
	manipulator X(angles);
	
	VectorXd x = TRANS(X.forwardKine(angles,3)); 
	VectorXd x2 = TRANS(X.forwardKine(angles,3));
	Vector3d dx;// std::vector< Vector3d >


	for(t=0.0;t<=TMAX;t+=twidth){
		angles << 2*M_PI, 2*M_PI, 2*M_PI*t/TMAX;		
			T=X.forwardKine(angles,3);
			x=TRANS(T);
			
		ofs << t << x.transpose() << "\n" << std::endl;
		dx=(x-x2)/twidth;			
		x2=x;
		
		dangles=(angles-angles2)/twidth;
		angles2=angles;
		
		J=X.getJacobian(angles);
		
		v=J*dangles;

// 		cout << "debug fk" << endl << X.forwardKine(Vector3d::Zero(),3) << endl ;	
		
// 		cout << "x" << endl << x << endl ;	
// 		
// 		cout << "x2" << endl << x2 << endl ;	
		
// 		cout << "dangles" << endl << dangles << endl ;	
		
		cout << "confirm Jacobian Effect ----> v,dx:hand velocity v:dx" <<endl<< v.head(3).transpose() <<endl<<  dx.transpose()  << endl;		//confirm Jacobian Effect
		
	}

/*	
		angles << 2*M_PI*0.001/TMAX, 2*M_PI, 2*M_PI;		
		T=X.forwardKine(angles,3);
		x=TRANS(T);	
		dx=(x-x2)/twidth;
*/	
	

/*	
	cout << "Jacobian\n" << J << endl;	
	
	J=X.testGetJacobian(dangles,dx);
	
	cout << "testJacobian\n" << J << endl;*/
	
	
	return 0;
}
