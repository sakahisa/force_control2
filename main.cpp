#include <iostream>
#include "Eigen/Dense"

using namespace Eigen;

Vector2d myFunction(Vector2d in)
{
	return in+in;
}


int main()
{
MatrixXd m(2,2);
m(0,0) = 3;
m(1,0) = 2.5;
m(0,1) = -1;
m(1,1) = m(1,0) + m(0,1);
VectorXd cx(2);
Vector2d c = Vector2d::Zero();
Vector2d d = Vector2d::Zero();
Matrix3d e = Matrix3d::Identity();
Matrix3d f = Matrix3d::Random();
Vector2d g = Vector2d::Zero();
g << 1.0, 1.0;
std::cout << ((m*m + m)*cx).transpose()<< std::endl;
std::cout << d<< std::endl;
std::cout << e.inverse()<< std::endl;
std::cout << f.inverse()<< std::endl;
std::cout << std::endl;
std::cout << g << std::endl;
std::cout << myFunction(g) << std::endl;
}
