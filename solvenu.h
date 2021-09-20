#ifndef SOLVENU_H
#define SOLVENU_H

#define PRINT_MAT(X) std::cout << #X << ":\n" << X << std::endl << std::endl
using namespace Mymath;

enum LIMITFlag{
    LIMITOFF=0,
    LIMITON=1,
    LIMITCHECK=2,
};

enum SolvFLAG{
  JACOBI=0,
  NEWTON=1,
  STEEPEST=2,
  DEFOKO=3,
};

enum JudgeFLAG{
    ZERO=0,
    DIFFZERO=1,
};

class Solvenu : public Funcvec{
    protected:
        long countlimit;
        VectorXd targetfx;
        VectorXd x;
        VectorXd upperlimit;
        VectorXd lowerlimit;
        LIMITFlag limitfl;
        double threshold;
    public:
        Solvenu();
        void setcountlimit(long a);
        void setlimit(VectorXd uplimit,VectorXd lowlimit);
        void setthreshold(double l_threshold);
        void unsetlimit();
        int checklimit(VectorXd &x);
        VectorXd sigmoidlimit(VectorXd &x,double alpha);
        VectorXd penaltyfunc(VectorXd &x);
        void settargetfx(VectorXd tfx);
        VectorXd gettargetfx();
        VectorXd function(VectorXd x) override;
        virtual VectorXd funcorg(VectorXd x);
        VectorXd functionerror(VectorXd x);
        VectorXd solve(VectorXd intx,SolvFLAG slflag,double l_alpha=0,JudgeFLAG jdgfl=ZERO);
        VectorXd newtonsolve(VectorXd intx,JudgeFLAG jdgfl);
        VectorXd newtonsolve(VectorXd intx,MatrixXd &l_jacobi,JudgeFLAG jdgfl);
        VectorXd steepsetdescentsolve(VectorXd intx,JudgeFLAG jdgfl);
        VectorXd steepsetdescentsolve(VectorXd intx,double l_alpha,JudgeFLAG jdgfl);
        ~Solvenu();
};

#endif