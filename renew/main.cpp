#include <iostream>
#include <cmath>
#include <fstream>
#include "Eigen/Dense"

#define twidth 0.001
#define TMAX 5.0
#define l1 1.0
#define l2 1.0

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
	Vector4d x;
	
	manipulator(){};
	manipulator(double a);
};
Link::Link(double l,double m,double n,double o){
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
	T(1,0)=sin(theta)*cos(alpha);		T(1,1)=cos(theta)*cos(alpha);	T(1,2)=-sin(alpha);		T(0,3)=-sin(alpha)*offset;
	T(2,0)=sin(theta)*sin(alpha);		T(2,1)=cos(theta)*sin(alpha);	T(2,2)=cos(alpha);		T(0,3)=cos(alpha)*offset;
	T(3,0)=0.0;							T(3,1)=0.0;						T(3,2)=0.0;				T(3,3)=1.0;
	return T;
}



manipulator::manipulator(double a){
	x(0)=a;
	x(1)=0.0;
	x(2)=0.0;
	x(3)=1.0;
}

/*
Matrix4d ForwardKine(double angle)
{
	Link ;
	Matrix4d T;
	Matrix4d F;
	
	T=a.Transform(angle);
	
	
	
	
	
	
	return F;
//todo
}
*/
int main()
{
	std::ofstream ofs("position.txt");
	
	double t;
	Matrix4d T_1,T_2;
	Link A(0.0,0.0,0.0,0.0);
	Link B(l1,0.0,0.0,0.0);
	A.a=l1;
	B.a=l2;
	
	manipulator X(l1);
	manipulator Y(l2);
	
	Vector4d m,n;
	
	for(t=0.0;t<TMAX;t+=twidth){
		A.theta=2*M_PI*t;
		B.theta=2*M_PI*t;
		T_1=A.Transform(A.theta);
		T_2=B.Transform(B.theta);
		m=T_1*X.x;
		n=m+T_1*T_2*Y.x;
		
		ofs << t << " " << m(0) << " " << m(1) << " " << m(2) << " " << n(0) << " " << n(1) << " " << n(2) << "\n" << std::endl;
	}
	return 0;
}
