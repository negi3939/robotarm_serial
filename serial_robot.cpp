#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <string>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/ioctl.h>
#include <fcntl.h>
#include <termios.h>
#include <unistd.h>
#include <cstdlib>
#include <iomanip>
#include <sys/time.h>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/StdVector>
#include <Eigen/LU>

#include "Negi39IKID.h"
#include "keyboard.h"
#include "serial.h"

class Robotserial{
    protected:
        Serial *ser;
        uint8_t *read8t;
        uint8_t *send8t;
        int motornum;
        int sendbyte;
        int readbyte;
        std::string devname;
        int boardrate;
        int *data;
    public:
        Robotserial();
        int setserial();
        int sendangle(VectorXd l_angle);
};

Robotserial::Robotserial(){
    motornum = 7;
    readbyte = 3;
    sendbyte= 1+2*motornum+1;
    read8t = new uint8_t[readbyte];
    send8t = new uint8_t[sendbyte];
    data = new int[motornum];
    setserial();
}

int Robotserial::setserial(){
    devname = "/dev/ttyS11";
    boardrate = B115200;
    
    ser = new Serial(boardrate,devname);
    send8t[0] = 0x0b;
    for(int ii=1;ii<sendbyte-1;ii++){
        send8t[ii] = 0x00;
    }
    ser->setdisplay();
}

int Robotserial::sendangle(VectorXd l_angle){
    if(l_angle.size()==(motornum-1)){
        data[motornum-1] = 70;
    }
    else if(l_angle.size()!=motornum){
        std::cout << "l_angle.size()\t" << l_angle.size() << std::endl;
        return 0;
    }
    
    for(int ii=0;ii<l_angle.size();ii++){
        data[ii] = (int)l_angle(ii);
        std::cout << std::dec << data[ii] << "," << std::flush;
    }
    std::cout << data[motornum-1] << std::endl;
    
    for(int ii=0;ii<motornum;ii++){
            send8t[2*ii+1] = ((uint16_t)data[ii] >> 8 ) & 0xFF; 
            send8t[2*ii+2] = (uint16_t)data[ii] & 0xFF;
    }
    ser->write_wcrc(send8t,sendbyte);
    return ser->read_get(read8t,readbyte,5);
}

int main(){
    keyboard ky;
    Robotserial rs;
    int jointnum = 6;
    struct timeval start_time,end_time;
    double time_second = 0.0;

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
    angle_deg_defo << 0.0d , -30.0d , -30.0d , -90.0d , -30.0d , 90.0d;
    gettimeofday(&start_time, NULL);
    while(1){
        rs.sendangle(angle_deg_defo);
        gettimeofday(&end_time, NULL);
        time_second = (double)(end_time.tv_sec - start_time.tv_sec)+(double)(end_time.tv_usec - start_time.tv_usec)*0.000001;
        if(time_second>1.0){break;}
    }
    roboik->setdefoko(angle_deg_defo,DEG);
    Vector3d pos;
    pos << -0.0,0.0,0.03; 
    double eul = -0.0d*M_PI/180.0d;
    ang = roboik->solve_relative(pos,YAW,eul,DEG);
    showvec(ang);
    sleep(1);
    gettimeofday(&start_time, NULL);
    while(1){
        rs.sendangle(ang);
        gettimeofday(&end_time, NULL);
        time_second = (double)(end_time.tv_sec - start_time.tv_sec)+(double)(end_time.tv_usec - start_time.tv_usec)*0.000001;
        if(time_second>1.0){break;}
    }
    delete roboik;
}