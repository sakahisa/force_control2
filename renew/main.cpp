#include <iostream>
#include <cmath>
#include <fstream>
#include "Eigen/Dense"
#include <vector>

#define twidth 0.001
#define TMAX 5.0
#define l1 1.0
#define l2 0.5
#define l3 0.25

#define TRANS(x) x.block<3,1>(0,3)
#define ROT(x) x.block<3,3>(0,0) 

using namespace Eigen;

class Link
{
	public:
	double a,alpha,offset,theta;
	
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
	MatrixXd testGetJacobian(VectorXd angle, Vector4d v);

	std::vector<Link* > links;
	
};


Link::Link(double l,double m,double n,double o)
{
	a=l;
	alpha=m;
	offset=n;
	theta=o;
}

Matrix4d Link::Transform(double angle)
{
	Matrix4d T;
	
	theta=angle;
	
	T(0,0)=cos(theta);					T(0,1)=-sin(theta);				T(0,2)=0.0;				T(0,3)=a;
	T(1,0)=sin(theta)*cos(alpha);		T(1,1)=cos(theta)*cos(alpha);	T(1,2)=-sin(alpha);		T(1,3)=-sin(alpha)*offset;
	T(2,0)=sin(theta)*sin(alpha);		T(2,1)=cos(theta)*sin(alpha);	T(2,2)=cos(alpha);		T(2,3)=cos(alpha)*offset;
	T(3,0)=0.0;							T(3,1)=0.0;						T(3,2)=0.0;				T(3,3)=1.0;

	return T;
}

manipulator::manipulator()
{
	int dof = 3;
	VectorXd l(3);
	l << 0.0, l1, l2;
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
	l << 0.0, l1, l2;
	VectorXd m(3);
	m << 0.0, 0.0, 0.0;
	VectorXd n(3);
	n << 0.0, 0.0, 0.0;
	VectorXd o(3);
	o = angle;
	
	for(int i = 0; i < dof; i++) links.push_back(new Link(l(i),m(i),n(i),o(i)));
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

		T *= links[i]->Transform(angle(i) ); // A->B = (*A).B
	}

	return T;
}

MatrixXd manipulator::getJacobian(VectorXd angle)
{
	MatrixXd Jacobian(6, angle.size());
	Vector3d z, p, p_minus1, pn;
	
	pn = TRANS(forwardKine(angle));
	
	for(int i = 0; links.size(); i++)
	{
		p_minus1 = TRANS(forwardKine(angle, i));
		p = pn - p_minus1;
		
		z = ROT(forwardKine(angle, i)).col(2);
		
		Jacobian.block<3,1>(0,i) = z.cross(p);
		Jacobian.block<3,1>(3,i) = z; 
	}
	
	return Jacobian;
}

MatrixXd manipulator::testGetJacobian(VectorXd angle, Vector4d v)
{
	MatrixXd Jacobian(6,angle.size());
	VectorXd dangle;
	dangle=angle;
	double V;
	
	for(int i; i<3; i++)
	{
		for(int j; j<3; j++)
		{
			Jacobian(i,j)=v(i)/angle(j);
		}
	}
	return Jacobian;
}

void initialize(Vector4d a[3])
{
	for(int i = 0; i < 3; i ++)
	{
		a[i](0)=0.0;
		a[i](1)=0.0;
		a[i](2)=0.0;
		a[i](3)=1.0;
	}
	a[0](0)=l1;
	a[1](0)=l2;
	a[2](0)=l3;
}
 
 
 
 
int main()
{
	std::ofstream ofs("position.txt");
	
	int i;
	double t;
	Link A;
	manipulator X;
	
	VectorXd angles(3);
	
	Vector4d a[3]=VectorXd::Zero(3);
	Vector4d a2[3]=VectorXd::Zero(3);
	Vector4d da[3]=VectorXd::Zero(3);
	
	initialize(a);
	
	
	
	for(t=0.0;t<TMAX;t+=twidth){
		angles << 2*M_PI*t/10, 2*M_PI*t/5, 2*M_PI*t;
		manipulator X(angles);
		for(i=0;i<3;i++)
		{
			a[i]=X.forwardKine(angles,i+1)*a[i];
			a[i](3)=1.0;
		}
		ofs << t << "\t" << a[0].transpose() << "\t" << a[1].transpose() << "\t" << a[2].transpose() << "\n" << std::endl;
		
		for(i=0;i<3;i++)
		{
			da[i]=(a[i]-a2[i])/twidth;			//verocity
			a2[i]=a[i];
		}
		initialize(a);
	}
	std::cout << X.getJacobian(angles) << " " << X.testGetJacobian(angles/t,da[2]) << "\n";
	return 0;
}
