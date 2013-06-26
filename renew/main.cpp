#include <iostream>
#include <cmath>
#include <fstream>
#include "Eigen/Dense"
#include <vector>

#define twidth 0.001
#define TMAX 5.0
#define l1 1.0
#define l2 1.0

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
	Vector4d x;
	
	manipulator(){};
	
	// camelCase	
	Matrix4d forwardKine(VectorXd angle);
	Matrix4d forwardKine(VectorXd angle, int to_idx);
	
	MatrixXd getJacobian(VectorXd angle);
	MatrixXd testGetJacobian(VectorXd angle);
	
	manipulator(double a);
	std::vector<Link* > links;
	
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
	//links.size() -> 0
	//links[0] = asdasd
	int dof = 3;
	double l[3] = {0.1, 0.05, 0.025};
	double m[3] = {0.1, 0.05, 0.025};
	double n[3] = {0.1, 0.05, 0.025};
	double o[3] = {0.1, 0.05, 0.025};
	
	for(int i = 0; i < dof; i++) links.push_back(new Link(l[i],m[i],n[i],o[i]));
}

Matrix4d manipulator::forwardKine(VectorXd angle)
{
	return forwardKine(angle, links.size());
	
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

MatrixXd manipulator::testGetJacobian(VectorXd angle)
{
	// TODO
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
	
	VectorXd angles(3);
	angles << 0, 0 , 0;
	
	std::cout << X.forwardKine(angles) << std::endl;
// 	forwardKine(angle, to_idx)
	return 0;
}
