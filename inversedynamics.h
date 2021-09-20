#ifndef INVDYN_H
#define INVDYN_H
#include "solvenu.h"
#include "inversekinematics.h"

#define PRINT_MAT(X) std::cout << #X << ":\n" << X << std::endl << std::endl
using namespace Mymath;

class invdSolvenu : public invkSolvenu {
  protected:
  public:
    invdSolvenu();
    invdSolvenu(int num);
    VectorXd funcorg(VectorXd x) override;
    VectorXd gettau(VectorXd f,VectorXd mom);
    void calcforce(VectorXd tau,Vector3d &f,Vector3d &mom);
    VectorXd getvel(VectorXd x,SolvFLAG solvfl);
    ~invdSolvenu();
};

#endif