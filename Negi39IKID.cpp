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

#include "Negi39IKID.h"


void DHparameter::set(double l_thetaoffset,double l_a,double l_d,double l_alpha){
    thetaoffset = l_thetaoffset;
    a = l_a;
    d = l_d;
    alpha = l_alpha;
}

double DHparameter::getthetaoffset(){return thetaoffset;}
double DHparameter::geta(){return a;}
double DHparameter::getd(){return d;}
double DHparameter::getalpha(){return alpha;}

Negi39IKID::Negi39IKID(){}

Negi39IKID::Negi39IKID(std::vector<DHparameter> l_dh){
    dh = l_dh;
    limitfl = LIMITOFF;
    init();
}

Negi39IKID::Negi39IKID(std::vector<DHparameter> l_dh,VectorXd l_uplimit,VectorXd l_lowlimit){
    dh = l_dh;
    uplimit = l_uplimit;
    lowlimit = l_lowlimit;
    limitfl = LIMITON;
    init();
}

void Negi39IKID::set(std::vector<DHparameter> l_dh){
    dh = l_dh;
    limitfl = LIMITOFF;
    init();
}

void Negi39IKID::set(std::vector<DHparameter> l_dh,VectorXd l_uplimit,VectorXd l_lowlimit){
    dh = l_dh;
    uplimit = l_uplimit;
    lowlimit = l_lowlimit;
    limitfl = LIMITON;
    init();
}

void Negi39IKID::init(){
    jointnum = dh.size();
    maninvk = new invkSolvenu(jointnum);
    maninvd = new invdSolvenu(jointnum);
    for(int ii=0;ii<dh.size();ii++){
        maninvk->setdhparameter(ii,dh[ii].getthetaoffset(),dh[ii].geta(), dh[ii].getd(),dh[ii].getalpha());
    }
    maninvd->copy(maninvk);
    //limit add
    if(limitfl!=LIMITOFF){
        maninvk->setlimit(uplimit,lowlimit);//可動範囲を設定（FLAGが立つ）
    }
    maninvk->setthreshold(0.00001d);//solverの判定thresholdを設定
    //init angle
    targetx = VectorXd::Zero(7);
    targetx(3) = 1.0d;
    angle_rad = VectorXd::Zero(6);
    angle_deg = VectorXd::Zero(6);
    angle_defo = VectorXd::Zero(6);
    //setdefoko();// defoko pose
    errorlimit = 0.001;
    maninvk->setcountlimit(1000);
    maninvk->setthreshold(sqrt((double)jointnum*0.001*0.001));

}

void Negi39IKID::setdefoko(VectorXd defo,ANGLEUNIT unit){
    switch (unit){
    case DEG:
        angle_deg = defo;
        angle_rad = degtorad(angle_deg);
        angle_defo = angle_rad;
        break;
    case RAD:
        angle_rad = defo;
        angle_deg = radtodeg(angle_rad);
        angle_defo = angle_rad;
        break;
    default:
        break;
    }
    setdefoko();
}

void Negi39IKID::setdefoko(){
    angle_rad = angle_defo;
    angle_deg = radtodeg(angle_defo);
    maninvk->calcaA(angle_rad,mattheta);
    //PRINT_MAT(mattheta);
    pos = mattheta.block(0,3,3,1);
    qua = maninvk->matrixtoquatanion(mattheta);
    targetx.block(0,0,3,1) = pos;
    targetx.block(3,0,4,1) = qua.block(0,0,4,1);
    std::cout << " target is  "<< std::flush;
    for(int ii=0;ii<6;ii++){
        std::cout << targetx(ii) << " , " << std::flush;
    }
    std::cout << targetx(6) << std::endl;
}

VectorXd Negi39IKID::degtorad(VectorXd deg){return M_PI/180.0d*deg;}
VectorXd Negi39IKID::radtodeg(VectorXd rad){return 180.0d/M_PI*rad;}

VectorXd Negi39IKID::solve_relative(Vector3d &target,EULERFL eulfl,double &eulangle,ANGLEUNIT unit){
    Matrix3d rotateR;
    rotateR = Matrix3d::Identity();
    for(int ii=0;ii<3;ii++){
        targetx(ii) += target(ii);     
    }
    switch (eulfl){
    case NON:
        break;    
    case ROLL:
        rotateR(1,1) = cos(eulangle);
        rotateR(1,2) = -sin(eulangle);
        rotateR(2,2) = cos(eulangle);
        rotateR(2,1) = sin(eulangle);
        break;
    case PITCH:
        rotateR(0,0) = cos(eulangle);
        rotateR(0,2) = sin(eulangle);
        rotateR(2,0) = -sin(eulangle);
        rotateR(2,2) = cos(eulangle);
        break;
    case YAW:
        rotateR(0,0) = cos(eulangle);
        rotateR(0,1) = -sin(eulangle);
        rotateR(1,1) = cos(eulangle);
        rotateR(1,0) = sin(eulangle);
        break;
    default:
        break;
    }

    if(eulfl!=NON){
        mattheta.block(0,0,3,3) = rotateR*mattheta.block(0,0,3,3);
    }
    pos = targetx.block(0,0,3,1);
    qua = maninvk->matrixtoquatanion(mattheta);
    targetx.block(3,0,4,1) = qua.block(0,0,4,1);
    //PRINT_MAT(mattheta);
    //targetx(1) = 0.0d;
    maninvk->settargetfx(targetx);
    angle_rad = maninvk->getangle(angle_rad,NEWTON);
    angle_deg = radtodeg(angle_rad);
    if(unit==DEG){return angle_deg;}
    return angle_rad; 
}

int Negi39IKID::solve_relative(Vector3d &target,EULERFL eulfl,double &eulangle,VectorXd &angle,ANGLEUNIT unit){//Solve 成功か否かを返す
    Matrix3d rotateR;
    int ret = 0;
    rotateR = Matrix3d::Identity();
    VectorXd old_targetx = targetx;
    for(int ii=0;ii<3;ii++){
        targetx(ii) += target(ii);     
    }
    switch (eulfl){
    case NON:
        break;    
    case ROLL:
        rotateR(1,1) = cos(eulangle);
        rotateR(1,2) = -sin(eulangle);
        rotateR(2,2) = cos(eulangle);
        rotateR(2,1) = sin(eulangle);
        break;
    case PITCH:
        rotateR(0,0) = cos(eulangle);
        rotateR(0,2) = sin(eulangle);
        rotateR(2,0) = -sin(eulangle);
        rotateR(2,2) = cos(eulangle);
        break;
    case YAW:
        rotateR(0,0) = cos(eulangle);
        rotateR(0,1) = -sin(eulangle);
        rotateR(1,1) = cos(eulangle);
        rotateR(1,0) = sin(eulangle);
        break;
    default:
        break;
    }

    if(eulfl!=NON){
        mattheta.block(0,0,3,3) = rotateR*mattheta.block(0,0,3,3);
    }
    pos = targetx.block(0,0,3,1);
    qua = maninvk->matrixtoquatanion(mattheta);
    targetx.block(3,0,4,1) = qua.block(0,0,4,1);
    //PRINT_MAT(mattheta);
    //targetx(1) = 0.0d;
    maninvk->settargetfx(targetx);
    angle_rad = maninvk->getangle(angle_rad,NEWTON);
    angle_deg = radtodeg(angle_rad);
    VectorXd ferr = maninvk->functionerror(angle_rad);
    for(int ii=0;ii<ferr.size();ii++){
        if(std::abs(ferr(ii))>errorlimit){
            ret+=1;
            break;
        }
    }
    if(ret>0){
        targetx = old_targetx;
        return ret;
    }else{
        if(unit==DEG){
            angle = angle_deg;
        }else{
            angle = angle_rad;
        }
    }
    return ret;
}

VectorXd Negi39IKID::solve_abstarg(VectorXd &target,ANGLEUNIT unit){
    for(int ii=0;ii<target.size();ii++){
        //if(ii!=1){
        targetx(ii) = target(ii);
        //}else{     
        //}
    }
    pos = targetx.block(0,0,3,1);
    qua.block(0,0,4,1) = targetx.block(3,0,4,1);
    //PRINT_MAT(mattheta);
    //targetx(1) = 0.0d;
    maninvk->settargetfx(targetx);
    angle_rad = maninvk->getangle(angle_rad,NEWTON);
    angle_deg = radtodeg(angle_rad);
    if(unit==DEG){return angle_deg;}
    return angle_rad;
}

int Negi39IKID::solve_abstarg(VectorXd &target,VectorXd &angle,ANGLEUNIT unit){//Solve 成功か否かを返す
    //maninvk->setcountlimit(100000);
    int ret = 0;
    for(int ii=0;ii<target.size();ii++){
        //if(ii!=1){
        targetx(ii) = target(ii);
        //}else{     
        //}
    }
    pos = targetx.block(0,0,3,1);
    qua.block(0,0,4,1) = targetx.block(3,0,4,1);
    //PRINT_MAT(mattheta);
    //targetx(1) = 0.0d;
    maninvk->settargetfx(targetx);
    angle_rad = maninvk->getangle(angle_rad,NEWTON);
    angle_deg = radtodeg(angle_rad);
    VectorXd ferr = maninvk->functionerror(angle_rad);
    for(int ii=0;ii<ferr.size();ii++){
        if(std::abs(ferr(ii))>errorlimit){
            ret+=1;
            break;
        }
    }
    //maninvk->setcountlimit(1000);
    if(ret>0){return ret;}
    if(unit==DEG){
        angle = angle_deg;
    }else{
        angle = angle_rad;
    }
    return ret;
}

int Negi39IKID::calc_foward(VectorXd &angle,VectorXd &posq,ANGLEUNIT unit){
    if(unit==DEG){
        angle_deg = angle;
    }else{
        angle_rad = angle;
        angle_deg = radtodeg(angle_rad);
    }
    angle_rad = degtorad(angle_deg);
    maninvk->calcaA(angle_rad,mattheta);
    pos = mattheta.block(0,3,3,1);
    qua = maninvk->matrixtoquatanion(mattheta);
    posq.block(0,0,3,1) = pos;
    posq.block(3,0,4,1) = qua.block(0,0,4,1);
}

void Negi39IKID::show_angle(ANGLEUNIT unit){
    VectorXd l_angle;
    switch (unit){
    case DEG:
        l_angle = angle_deg;
        break;
    case RAD:
        l_angle = angle_rad;
        break;
    default:
        break;
    }
    std::cout << "angles are \t";
    for(int ii=0;ii<maninvk->getjointnum()-1;ii++){
        std::cout << l_angle(ii) << " , ";
    }
    std::cout << l_angle(maninvk->getjointnum()-1) <<  std::endl;
}

VectorXd Negi39IKID::setangle(VectorXd ang){
    angle_deg = ang;
    angle_rad = degtorad(angle_deg);
    maninvk->calcaA(angle_rad,mattheta);
    //log_fs.write_fn(memo);PRINT_MAT(mattheta);
    pos = mattheta.block(0,3,3,1);
    qua = maninvk->matrixtoquatanion(mattheta);
    targetx.block(0,0,3,1) = pos;
    targetx.block(3,0,4,1) = qua.block(0,0,4,1);
    return targetx;
}

#if defined(NEGI_IS_MAIN)
int main(){
    int jointnum = 6;
    std::vector<DHparameter> mytobot(jointnum);
    mytobot[0].set( 0.0d*M_PI,     0.0d, 0.0965d, -0.5d*M_PI); 
    mytobot[1].set( -0.5d*M_PI, 0.1075d,    0.0d,  0.0d*M_PI);
    mytobot[2].set( 0.0d*M_PI,     0.0d,    0.0d,  0.5d*M_PI);
    mytobot[3].set( 0.0d*M_PI, 0.134d, 0.0d,  -0.5d*M_PI);
    mytobot[4].set( 0.5d*M_PI,   0.0d, -0.015d, 0.5d*M_PI);
    mytobot[5].set( 0.5d*M_PI, 0.015d, 0.205d, 0.0d*M_PI);
    
    Negi39IKID *roboik = new Negi39IKID(mytobot);
    VectorXd ang(jointnum);
    VectorXd angle_deg_defo(6);
    angle_deg_defo << 0.0d , -30.0d , -30.0d , -90.0d , -30.0d , 0.0d;
    roboik->setdefoko(angle_deg_defo,DEG);
    Vector3d pos;
    pos << -0.0,0.0,-0.0; 
    double eul = -0.0d*M_PI/180.0d;
    roboik->solve_relative(pos,YAW,eul);
    roboik->show_angle();
    delete roboik;
}
#endif