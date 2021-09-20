#ifndef INVKIN_H
#define INVKIN_H
#include "solvenu.h"

#define PRINT_MAT(X) std::cout << #X << ":\n" << X << std::endl << std::endl
using namespace Mymath;

class invkSolvenu : public Solvenu {
  protected:
    int jointnum;
    double *time;
    double *pre_time;
    Matrix4d *aA;
    Matrix4d *aAdis;
    Matrix4d *aAtheta;
    Matrix4d *aAthetaoff;
    Matrix4d *aAaal;
    Matrix4d *aAalp;
    double *aal;
    double *alp;
    double *dis;
    double *thetaoff;
    Matrix4d *aTt;
    Vector3d *rra;
    Vector3d *ppa;
    Vector3d *zz;
    MatrixXd *jacobi;
    void init();
  public:
    invkSolvenu();
    invkSolvenu(int num);
    void copy(const invkSolvenu &invk);
    void copy(invkSolvenu *invk);
    void setjointnum(int n);
    void settime(const double &t);
    void settime(double *t);
    void calcaTt();
    void calcjacobi();
    MatrixXd getjacobi();
    MatrixXd getpseudoinvjacobi();
    int getjointnum();
    double gettime();
    double* gettimead();
    double getaal(int n);
    double getalp(int n);
    double getdis(int n);
    double getthetaoff(int n);
    VectorXd funcorg(VectorXd x) override;
    void calcaA(VectorXd x);
    void calcaA(VectorXd x,Matrix4d &reta);
    void setdhparameter(int num,double thoff,double aa,double di,double alph);
    Vector4d matrixtoquatanion(Matrix4d mat);
    Vector4d matrixtoquatanion(Matrix3d mat);
    Matrix3d quataniontomatrix(Vector4d qua);
    VectorXd getangle(VectorXd x,SolvFLAG solvfl=DEFOKO);
    ~invkSolvenu();
};

#endif