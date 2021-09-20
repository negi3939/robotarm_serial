#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <unistd.h>
#include <getopt.h>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <stdint.h>
#include <vector>
#include <math.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/times.h>
#include <sys/time.h>
#include <sys/ioctl.h>
#include <fcntl.h>
#include <bitset>
#include <termios.h>
#include <pthread.h>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/StdVector>
#include <Eigen/LU>
#include "mymath.h"
#include "solvenu.h"

Solvenu::Solvenu(){countlimit=1000;limitfl=LIMITOFF;threshold=0.000001d;}
Solvenu::~Solvenu(){}

void Solvenu::setcountlimit(long a){countlimit = a;}
void Solvenu::setthreshold(double l_threshold){threshold = l_threshold;}

void Solvenu::setlimit(VectorXd uplimit,VectorXd lowlimit){
    upperlimit = uplimit;
    lowerlimit = lowlimit;
    limitfl = LIMITON;
}

void Solvenu::unsetlimit(){
    limitfl = LIMITOFF;
}

int Solvenu::checklimit(VectorXd &x){
    if(limitfl==LIMITOFF){return 0;}
    int fl=0;
    for(int ii=0;ii<x.size();ii++){
        //errval(ii) = step(x(ii) - upperlimit(ii))*abs(x(ii) - upperlimit(ii)) - step( - x(ii) + lowerlimit(ii))*abs( - x(ii) + lowerlimit(ii)); 
        if(x(ii)> upperlimit(ii)){x(ii) = upperlimit(ii);fl++;}
        else if(x(ii)< lowerlimit(ii)){x(ii) = lowerlimit(ii);fl++;}
    }
    return fl;
}

VectorXd Solvenu::sigmoidlimit(VectorXd &x,double alpha){
    VectorXd ans = VectorXd::Zero(1);
    if(limitfl==LIMITOFF){return ans;}
    for(int ii=0;ii<x.size();ii++){
        ans(0) += sigmoid(x(ii)-upperlimit(ii),alpha);
    }
    return ans;
}

VectorXd Solvenu::penaltyfunc(VectorXd &x){
    VectorXd ans =VectorXd::Zero(x.size()) ;
    for(int ii=0;ii<x.size();ii++){
        ans(ii) += step(x(ii)-upperlimit(ii))*(x(ii)-upperlimit(ii))*(x(ii)-upperlimit(ii));
        ans(ii) += step(-x(ii)+lowerlimit(ii))*(-x(ii)+lowerlimit(ii))*(-x(ii)+lowerlimit(ii));
    }
    return ans;
}

void Solvenu::settargetfx(VectorXd tfx){targetfx = tfx;}
VectorXd Solvenu::gettargetfx(){return targetfx;}

VectorXd Solvenu::function(VectorXd x){
    return functionerror(x);
}

VectorXd Solvenu::funcorg(VectorXd x){
    VectorXd ans(2);
    MatrixXd aA(ans.size(),x.size());
    aA << 1,0,0,1;
    ans = aA*x;
    return ans;
}

VectorXd Solvenu::functionerror(VectorXd x){
    //checklimit(x);
    VectorXd buf = funcorg(x);
    VectorXd trgfx(buf.size());
    if(limitfl!=LIMITON){ 
        return  funcorg(x) - targetfx;
    }else{
        //PRINT_MAT(funcorg(x));
        //PRINT_MAT(targetfx);
        //PRINT_MAT(funcorg(x)- targetfx);
        MatrixXd oone = 100.0d*MatrixXd::Ones(targetfx.size(),x.size());
        return  funcorg(x) - targetfx + oone*penaltyfunc(x);
    }
}

VectorXd Solvenu::solve(VectorXd intx,SolvFLAG slflag=NEWTON,double l_alpha=0,JudgeFLAG jdgfl=ZERO){
    switch (slflag){
    case NEWTON:
        return newtonsolve(intx,jdgfl);
        break;
    case STEEPEST:
        if(l_alpha<=0){return steepsetdescentsolve(intx,jdgfl);}
        return steepsetdescentsolve(intx,l_alpha,jdgfl);
        break;
    default:
        return intx;
        break;
    }
}

VectorXd Solvenu::newtonsolve(VectorXd intx,JudgeFLAG jdgfl){
    x = intx;
    long count = 0;
    VectorXd dx = intx;
    MatrixXd buf;
    VectorXd chk;
    while(1){
        dx = x;
        JacobiSVD<MatrixXd> svd(diffvec(x,this), ComputeThinU|ComputeThinV);
        x = x - svd.solve(functionerror(x));//limit ここを書き換える必要がある
        //PRINT_MAT(x);
        chk = MatrixXd::Ones(1,x.size())*absmat(x - dx);// + sigmoidlimit(x,1000);
        if(count>countlimit){
            //std::cout <<"CAUTIONCAUTIONCAUTIONCAUTIONCAUTION step is more than" << countlimit << "CAUTIONCAUTIONCAUTIONCAUTIONCAUTIONC"<<std::endl;
            //PRINT_MAT(functionerror(x));
            break;
        }
        if ((jdgfl==ZERO)&&(functionerror(x).norm() < threshold)){ break;}
        else if((jdgfl==DIFFZERO)&&(diffnorm(x,this).norm()<threshold)){break;}
        count++;
    }
    int checkret=0;
    if(limitfl==LIMITCHECK){checkret = checklimit(x);}
    if(checkret){std::cout << "solve fail :" << checkret << std::endl;}
    return x;
}

VectorXd Solvenu::newtonsolve(VectorXd intx,MatrixXd &l_jacobi,JudgeFLAG jdgfl){
    x = intx;
    long count = 0;
    VectorXd dx = intx;
    MatrixXd buf;
    VectorXd chk;
    while(1){
        dx = x;
        JacobiSVD<MatrixXd> svd(l_jacobi, ComputeThinU|ComputeThinV);
        x = x - svd.solve(functionerror(x));//limit ここを書き換える必要がある
        //PRINT_MAT(x);
        chk = MatrixXd::Ones(1,x.size())*absmat(x - dx);// + sigmoidlimit(x,1000);
        if(count>countlimit){
            //std::cout <<"CAUTIONCAUTIONCAUTIONCAUTIONCAUTION step is more than" << countlimit << "CAUTIONCAUTIONCAUTIONCAUTIONCAUTIONC"<<std::endl;
            //PRINT_MAT(functionerror(x));
            break;
        }
        if ((jdgfl==ZERO)&&(functionerror(x).norm() < threshold)){ break;}
        else if((jdgfl==DIFFZERO)&&(diffnorm(x,this).norm()<threshold)){break;}
        count++;
    }
    int checkret=0;
    if(limitfl==LIMITCHECK){checkret = checklimit(x);}
    if(checkret){std::cout << "solve fail :" << checkret << std::endl;}
    return x;
}

VectorXd Solvenu::steepsetdescentsolve(VectorXd intx,JudgeFLAG jdgfl){
    long count = 0;
    VectorXd x = intx;
    MatrixXd diffv;
    double gold_r = 0.6180339887d;
    double bottom_alpha,top_alpha,alpha1,alpha2;
    MatrixXd s;
    while(1){
        diffv =  diffnorm(x,this);
        if((jdgfl==ZERO)&&(functionerror(x).norm()<threshold)){break;}
        else if((jdgfl==DIFFZERO)&&(diffnorm(x,this).norm()<threshold)){break;}
        //std::cout << "norm is " <<functionerror(x).norm() << std::endl;
        s = -diffv.transpose();
        bottom_alpha=0.0d;
        top_alpha=10.0d;
        alpha1 = bottom_alpha + (1.0d -gold_r)*(top_alpha - bottom_alpha);
        alpha2 = bottom_alpha + gold_r*(top_alpha - bottom_alpha);
        while(1){
            if(functionerror(x + alpha1*s).norm()<functionerror(x + alpha2*s).norm()){
                top_alpha = alpha2;
                alpha2 = alpha1;
                alpha1 = bottom_alpha + (1.0d -gold_r)*(top_alpha - bottom_alpha);
            }else{
                bottom_alpha = alpha1;
                alpha1 = alpha2;
                alpha2 = bottom_alpha + gold_r*(top_alpha - bottom_alpha);
            }
            if(count>countlimit){
                //std::cout <<"CAUTIONCAUTIONCAUTIONCAUTIONCAUTION step is more than" << countlimit << "CAUTIONCAUTIONCAUTIONCAUTIONCAUTIONC"<<std::endl;
                //PRINT_MAT(functionerror(x));
                break;
            }
            //std::cout << "bottom_alpha : "<< bottom_alpha << " top_alpha : "<< top_alpha << std::endl;
            if(std::abs(bottom_alpha - top_alpha) <threshold){break;}
        }
        x = x + 0.5d*(bottom_alpha+top_alpha)*s;
    }
    int checkret=0;
    if(limitfl==LIMITCHECK){checkret = checklimit(x);}
    if(checkret){std::cout << "solve fail :" << checkret << std::endl;}
    return x;
}

VectorXd Solvenu::steepsetdescentsolve(VectorXd intx,double l_alpha,JudgeFLAG jdgfl){
    long count = 0;
    VectorXd x = intx;
    MatrixXd diffv;
    MatrixXd s;
    while(1){
        //diffv =  diffnorm(x,this);
        if((jdgfl==ZERO)&&(functionerror(x).norm()<threshold)){break;}
        else if((jdgfl==DIFFZERO)&&(diffnorm(x,this).norm() < threshold)){break;}

        //std::cout << "alppha is const. norm is " <<functionerror(x).norm() << std::endl;
        s = -diffv.transpose();
        x = x + l_alpha*s;
    }
    int checkret=0;
    if(limitfl==LIMITCHECK){checkret = checklimit(x);}
    if(checkret){std::cout << "solve fail :" << checkret << std::endl;}
    return x;
}