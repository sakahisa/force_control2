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

using namespace std;
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
	MatrixXd testGetJacobian(VectorXd angle, Vector3d v);

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
	int dof = 4;
	VectorXd l(4);
	l << 0.0, l1, l2, l3;
	VectorXd m(4);
	m << 0.0, 0.0, 0.0, 0.0;
	VectorXd n(4);
	n << 0.0, 0.0, 0.0, 0.0;
	VectorXd o(4);
	o << 0.0, 0.0, 0.0, 0.0;
	
	for(int i = 0; i < dof; i++) links.push_back(new Link(l(i),m(i),n(i),o(i)));
}

manipulator::manipulator(VectorXd angle)
{
	int dof = 4;
	VectorXd l(4);
	l << 0.0, l1, l2, l3;
	VectorXd m(4);
	m << 0.0, 0.0, 0.0, 0.0;
	VectorXd n(4);
	n << 0.0, 0.0, 0.0, 0.0;
	VectorXd o(4);
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

MatrixXd manipulator::testGetJacobian(VectorXd angle, Vector3d v)
{
	MatrixXd Jacobian(6,angle.size());
	
	for(int i; i<3; i++)
	{
		for(int j; j<angle.size(); j++)
		{
			Jacobian(i,j) = v(i) / angle(j);
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
	manipulator X;
	
	VectorXd angles(4);
	VectorXd angles2(4);
	VectorXd dangles(4);
	
	Vector3d a[4]; 
	Vector3d a2[4];
	Vector3d da[4];
	
	VectorXd omega;
	VectorXd v(6);
	MatrixXd J(6,angles.size());
	Matrix4d T;
	
	for(t=0.0;t<TMAX;t+=twidth){
		angles << 2*M_PI*t/TMAX, 2*M_PI, 2*M_PI, 0.0;
		manipulator X(angles);
		for(i=0;i<3;i++)
		{
			T=X.forwardKine(angles,i+2);
			a[i]=TRANS(T);
		}
		ofs << t << "\t" << a[0].transpose() << "\t" << a[1].transpose() << "\t" << a[2].transpose() << "\n" << std::endl;
		for(i=0;i<3;i++)
		{
			if(t>twidth/2)
			{
				da[i]=(a[i]-a2[i])/twidth;			
			}
			a2[i]=a[i];
		}
		
		if(t>twidth/2)
		{
			dangles=(angles-angles2)/twidth;
		}
		angles2=angles;
		
	}
	
	J=X.getJacobian(angles);
	
	v=J*dangles;
	
	cout << "confirm Jacobian Effect ----> v,da[2]:hand velocity \nv:\n" << v << "\nda[2]:\n" << da[2]  << endl;		//confirm Jacobian Effect
	
	cout << "Jacobian\n" << J << endl;	
	
	J=X.testGetJacobian(dangles,da[2]);
	
	cout << "testJacobian\n" << J << endl;
	
	
	return 0;
}
