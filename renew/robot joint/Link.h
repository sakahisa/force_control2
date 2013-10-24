#ifndef LINK_H
#define LINK_H

#include "Eigen/Dense"

using namespace Eigen;

class Link
{
	public:
	double a, alpha, offset, theta_offset, mass;
	
	Link(){};
	Link(double l, double m, double n, double o, double p);
	Matrix4d Transform(double angles);
};

#endif //LINK_H
