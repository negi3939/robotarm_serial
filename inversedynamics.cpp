#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <unistd.h>
#include <getopt.h>
#include <stdlib.h>
#include <string>
#include <stdio.h>
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
#include "inversekinematics.h"
#include "inversedynamics.h"

invdSolvenu::invdSolvenu():invkSolvenu(){}
invdSolvenu::invdSolvenu(int num):invkSolvenu(num){}
invdSolvenu::~invdSolvenu(){delete jacobi;delete[] ppa;delete[] rra;delete[] aTt;delete[] zz;}

VectorXd invdSolvenu::funcorg(VectorXd x){
    calcaTt();
    calcjacobi();
    return (*jacobi)*x;
}

VectorXd invdSolvenu::gettau(VectorXd f,VectorXd mom){
    calcaTt();
    calcjacobi();
    VectorXd fm,tau;
    fm.resize(f.size()+mom.size(),1);
    fm.block(0,0,f.size(),1) = f;
    fm.block(f.size(),0,mom.size(),1) = mom;
    tau = jacobi->transpose()*fm;
    return tau;
}

void invdSolvenu::calcforce(VectorXd tau,Vector3d &f,Vector3d &mom){
    calcaTt();
    calcjacobi();
    VectorXd fmom(6);
    fmom.block(0,0,3,1) = f;
    fmom.block(3,0,3,1) = mom;
    /*JacobiSVD<MatrixXd> svd(jacobi->transpose(),ComputeThinU|ComputeThinV);
    fmom = svd.solve(tau);*/
    fmom = pseudo_inv(jacobi->transpose())*tau;
    f = fmom.block(0,0,3,1);
    mom = fmom.block(3,0,3,1);//動作未確認*/
}

VectorXd invdSolvenu::getvel(VectorXd x,SolvFLAG solvfl=DEFOKO){
    VectorXd vel(x.size());
    switch (solvfl){
        case NEWTON://ニュートン法で解く（遅い）
            vel = solve(x,NEWTON);
            break;
        case STEEPEST://最急降下法で解く（速い）
            vel = solve(x,STEEPEST);
            break;
        case DEFOKO:
            if((jointnum>6)||(limitfl==LIMITON)){//もし関節数が7以上もしくは関節可動範囲制約がある
                vel = solve(x,STEEPEST);//最急降下法で解く
            }else{
                vel = solve(x,NEWTON);//ニュートン法で解く
            }
        default:
            break;
    }
    return vel;
}

#if defined(ID_IS_MAIN)
VectorXd forward_dynamics(invdSolvenu *maninvd,VectorXd &angle,Vector3d &forcev,Vector3d &momentv){
    VectorXd tau;
    angle(0) = -0.25d*M_PI;
    angle(1) = -0.25d*M_PI;
    angle(2) = -0.25d*M_PI;
    angle(3) = -0.25d*M_PI;
    angle(4) = -0.25d*M_PI;
    angle(5) = -0.25d*M_PI;
    angle(6) = -0.25d*M_PI;
    maninvd->calcaA(angle);
    forcev(0) = 1.0d;
    forcev(1) = 1.0d;
    forcev(2) = 1.0d;
    momentv(0) = 0.0d;
    momentv(1) = 0.0d;
    momentv(2) = 0.0d;
    tau = maninvd->gettau(forcev,momentv);
    PRINT_MAT(tau);
    return tau;
}

void inverse_dynamics(invdSolvenu *maninvd,VectorXd &angle,VectorXd &ctauv,Vector3d &forcev,Vector3d &momentv){
    double torqconstant = 3.70d;
    angle(0) = -0.25d*M_PI;
    angle(1) = -0.25d*M_PI;
    angle(2) = -0.25d*M_PI;
    angle(3) = -0.25d*M_PI;
    angle(4) = -0.25d*M_PI;
    angle(5) = -0.25d*M_PI;
    angle(6) = -0.25d*M_PI;
    maninvd->calcaA(angle);
    ctauv(0) = 1.0d;
    ctauv(1) = 1.0d;
    ctauv(2) = 1.0d;
    ctauv(3) = 1.0d;
    ctauv(4) = 1.0d;
    ctauv(5) = 1.0d;
    ctauv(6) = 1.0d;
    maninvd->calcforce(torqconstant*ctauv,forcev,momentv);
    PRINT_MAT(forcev);
    PRINT_MAT(momentv);
}

void check_jacobi(invkSolvenu *maninvk,invdSolvenu *maninvd,VectorXd angle,Matrix4d mattheta){
    Vector4d qua;
    Vector3d pos,old_pos;
    VectorXd vel(6);
    VectorXd angular_vel(7);
    VectorXd jacobi_vel(7);
    int old_pos_fl= 0;
    VectorXd targetx(7);
    VectorXd old_angle = angle;
    double time = 0.0,deltat = 0.1;
    angle(0) = -0.25d*M_PI;
    angle(1) = -0.25d*M_PI;
    angle(2) = -0.25d*M_PI;
    angle(3) = -0.25d*M_PI;
    angle(4) = -0.25d*M_PI;
    angle(5) = -0.25d*M_PI;
    angle(6) = -0.25d*M_PI;
    old_angle = angle;
    maninvk->calcaA(angle,mattheta);
    pos = mattheta.block(0,3,3,1);
    old_pos = pos;
    while(pos(1)>0.44){
        pos(1) -= 0.001;
        for(int ii=0;ii<3;ii++){
            vel(ii) = (pos(ii) - old_pos(ii))/deltat;
        }
        for(int ii=3;ii<6;ii++){
            vel(ii) = 0;
        }
        qua = maninvk->matrixtoquatanion(mattheta);//回転行列からクオータニオンへ変換
        std::cout << " x : " << pos(0) << " y : " << pos(1) << " z : " << pos(2) << std::endl;
        targetx.block(0,0,3,1) = pos;
        targetx.block(3,0,4,1) = qua.block(0,0,4,1);
        maninvk->settargetfx(targetx);
        angle = maninvk->getangle(angle,DEFOKO);
        maninvd->calcaA(angle);
        maninvd->settargetfx(vel);
        std::cout << "angle_vels are \t";
        for(int ii=0;ii<maninvk->getjointnum()-1;ii++){
            angular_vel(ii) = (angle(ii) - old_angle(ii))/deltat;
            std::cout << angle(ii) << " , ";
        }
        std::cout << angle(maninvk->getjointnum()-1) <<  std::endl;
        //PRINT_MAT(maninvd->getpseudoinvjacobi()*maninvd->getjacobi());
        jacobi_vel.block(0,0,6,1) = maninvd->getjacobi()*angular_vel;//maninvd->getvel(jacobi_vel);
        std::cout << "jacobi_vels are \t";
        for(int ii=0;ii<maninvd->getjointnum()-1;ii++){
            std::cout << vel(ii) - jacobi_vel(ii) << " , ";
        }
        std::cout << jacobi_vel(maninvk->getjointnum()-1) <<  std::endl;
        old_pos = pos;
        old_angle = angle;
        time += deltat;
        old_pos_fl = 1;
    }
}

int main(){
    int ii,jointn = 7;
    invdSolvenu *maninvd;
    invkSolvenu *maninvk;
    maninvd = new invdSolvenu(jointn);
    maninvk = new invkSolvenu(jointn);
    /*RT CRANE*/
    maninvk->setdhparameter(0,0.0d*M_PI,0.0d,0.064d,-0.5d*M_PI);//(int num,double thoff,double aa,double di,double alph);
    maninvk->setdhparameter(1,0.0d*M_PI,0.0d,0.0d,0.5d*M_PI);//(int num,double thoff,double aa,double di,double alph);
    maninvk->setdhparameter(2,0.0d*M_PI,0.0d,0.065d+0.185d,-0.5d*M_PI);//(int num,double thoff,double aa,double di,double alph);
    maninvk->setdhparameter(3,0.0d*M_PI,0.0d,0.0d,0.5d*M_PI);//(int num,double thoff,double aa,double di,double alph);
    maninvk->setdhparameter(4,0.0d*M_PI,0.0d,0.121d+0.129d,-0.5d*M_PI);//(int num,double thoff,double aa,double di,double alph);
    maninvk->setdhparameter(5,0.0d*M_PI,0.0d,0.0d,0.5d*M_PI);//(int num,double thoff,double aa,double di,double alph);
    maninvk->setdhparameter(6,0.0d*M_PI,0.0d,0.019d+0.084d,0.0d);//(int num,double thoff,double aa,double di,double alph);
    maninvd->copy(maninvk);
    //limit add
    VectorXd uplimit(7);
    VectorXd lowlimit(7);
    uplimit <<    2.72271 ,  0.5*M_PI  ,   2.72271  ,        0 ,   2.72271 ,    0.5*M_PI ,   2.89725;//可動上限範囲を設定
    lowlimit <<  -2.72271 , -0.5*M_PI  ,  -2.72271  , -2.79253 ,  -2.72271 ,   -0.5*M_PI ,  -2.89725;//可動下限範囲を設定
    maninvk->setlimit(uplimit,lowlimit);//可動範囲を設定（FLAGが立つ）
    VectorXd angle = VectorXd::Zero(jointn);//joint angle
    VectorXd ctauv = VectorXd::Zero(jointn);//current
    Matrix4d mattheta = MatrixXd::Identity(4,4);//回転変位行列
    Vector3d forcev,momentv;//手先力,手先モーメント

    /*test*/
    check_jacobi(maninvk,maninvd,angle,mattheta);
    //forward_dynamics(maninvd,angle,forcev,momentv);//calc FD test
    //inverse_dynamics(maninvd,angle,ctauv,forcev,momentv);//calc ID test

    delete maninvd;
    delete maninvk;
    return 0;
}
#endif