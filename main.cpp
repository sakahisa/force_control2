#include <iostream>
#include <cmath>
#include "Eigen/Dense"

#define twidth 0.001
#define TMAX 5.0

using namespace Eigen;

class link
{
	public:
	double a, alpha, offset, theta;

	link(){};
	Matrix4d Transform(double angle);
};

class manipulator
{
	private:
	Vector4d x;
	Vector4d v;
	Vector4d a;
	manipulator();
	manipulator(Vector4d ac);
};

Matrix4d link::Transform(double angle)
{
	Matrix4d T;
	int i,j;
	
	T(0,3)=0.0;
	T(1,3)=a*cos(angle);
	T(2,3)=a*sin(angle);
	T(3,3)=1.0;
	for(i=0;i<3;i++){
		T(3,i)=0.0;
	}
	T(0,0)=1.0;		T(0,1)=0.0;			T(0,2)=0.0;
	T(1,0)=0.0;		T(1,1)=cos(angle);	T(1,2)=-sin(angle);
	T(2,0)=0.0;		T(2,1)=sin(angle);	T(2,2)=cos(angle);
	return T;
};


Matrix4d ForwardKine(double angle)
{
	Matrix4d F;
	
	
	
	
	
	return F;
//todo
}
void manipulator::manipulator(){
	a=Zero(4);
	v=Zero(4);
}

void manipulator::manipulator(Vector4d ac){
	a=ac;
	v+=ac*twidth;
	x+=v*twidth;
}


int main()
{
	double t;
	double the_1,the_2;
	Matrix4d F_1,F_2;
	link A,B;
	manipulator X,Y;
	manipulator cmd;
	
	for(t=0.0;t<TMAX;t+=twidth)
		F_1=Forward(the_1);
		F_2=Forward(the_2);
		
		
		
}
