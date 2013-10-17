#include "Eigen/Dense"
#include <vector>
#include <iostream>

using namespace std;
using namespace Eigen;

typedef Matrix<double, 6, 1> Vector6d;

class baseClass
{
    public:
    Vector6d Dis;
    void LPF(Vector6d forceDis);
};

class Link
{
    public:
    double a, alpha, offset, theta_offset, mass;

    Link(){};
    Link(double l,double m,double n,double o, double p);
    Matrix4d Transform(double angle);
};

class force : public baseClass
{
    public:
    Vector6d Cmd,Res,Dis,dDis;

    force();
};

class tau
{
    public:
    VectorXd Cmd,Res,Dis;

    tau();
    void Disturbance(VectorXd omega);
};


class manipulator : public baseClass
{
    public:

    force F;
    tau T;

    manipulator();

    Matrix4d forwardKine(VectorXd angle);
    Matrix4d forwardKine(VectorXd angle, int to_idx);

    MatrixXd getJacobian(VectorXd angle);

    Vector6d torque2Force(VectorXd angles, VectorXd tauRef);
    VectorXd force2Torque(VectorXd angles, Vector6d forceRef);
    VectorXd forceControl(VectorXd angles, Vector6d forceRef, VectorXd omega);

    MatrixXd pseudoInverse(MatrixXd Jacobian);
    VectorXd moment(VectorXd angles, Vector3d fRef);

    VectorXd invKine(Vector3d Error, VectorXd angles);

    Vector4d COMpos(int start, VectorXd angles);
    MatrixXd getCOMJacobian(VectorXd angles);

    std::vector<Link* > links;
};
