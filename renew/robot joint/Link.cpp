#include "Link.h"

Link::Link(double l,double m,double n,double o, double p)
{
	a = l;
	alpha = m;
	offset = n;
	theta_offset = o;
	mass = p;
}
//Transform Matrix
Matrix4d Link::Transform(double angles)
{
	Matrix4d T;
	
	double theta = angles + theta_offset;
	
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
