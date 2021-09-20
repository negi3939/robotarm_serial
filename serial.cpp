#include <iostream>
#include <sstream>
#include <fstream>
#include <math.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/ioctl.h>
#include <fcntl.h>
#include <termios.h>
#include <unistd.h>
#include <cstdlib>
#include <iomanip>
#include <sys/time.h>


#include "keyboard.h"
#include "serial.h"


Serial::Serial(){
    crc8ccit_poly = 0x6F;
    baudRate = B115200;
    buf = new unsigned char[255];
    portname = new char[12];
    portname = "/dev/ttyS1";
    init();
}

Serial::Serial(int baundrate){
    crc8ccit_poly = 0x6F;
    baudRate = baundrate;
    buf = new unsigned char[255];
    init();
}
Serial::Serial(int baundrate,char *devname){
    crc8ccit_poly = 0x6F;
    buf = new unsigned char[255];
    int sizename = 30;
    portname = new char[sizename];
    for(int ii=0;ii<sizename;ii++){
        portname[ii] = devname[ii];
    }
    std::cout << portname << std::endl;
    baudRate = baundrate;
    init();
}
Serial::Serial(int baundrate,std::string devname){
    crc8ccit_poly = 0x6F;
    buf = new unsigned char[255];
    int sizename = 30;
    portname = devname.c_str();
    std::cout << portname << std::endl;
    baudRate = baundrate;
    init();
}

Serial::Serial(int baundrate,std::string devname,uint8_t l_poly){
    crc8ccit_poly = l_poly;
    buf = new unsigned char[255];
    int sizename = 30;
    portname = devname.c_str();
    std::cout << portname << std::endl;
    baudRate = baundrate;
    init();
}

int Serial::init(){

    fd = open(portname, O_RDWR |  O_NONBLOCK);
    if (fd < 0){
        std::cout<< "open error!" <<std::endl;
        return -1;
    }
    tcgetattr(fd, &oldtio);
    bzero(&tio,sizeof(tio));
    
    tio.c_cflag |= CREAD;               // 受信有効
    tio.c_cflag |= CLOCAL;              // ローカルライン（モデム制御なし）
    tio.c_cflag |= CS8;                 // データビット:8bit
    tio.c_cflag |= 0;                   // ストップビット:1bit
    tio.c_cflag |= 0;                   // パリティ:None
    tio.c_oflag = 0;
    tio.c_lflag &= ~(ICANON | ECHO | ECHOE | ISIG);
  
    cfsetospeed( &tio, baudRate );
    cfsetispeed( &tio, baudRate );
    cfmakeraw(&tio);                    // RAWモード
    tcflush(fd,TCOFLUSH);
    tcflush(fd,TCIFLUSH);
    tcsetattr( fd, TCSANOW, &tio );     // デバイスに設定を行う
    ioctl(fd, TCSETA, &tio);

    gen_crc8ccit_table();
    display = SerialD::OFF;
    return 0;
}

int Serial::setdisplay(){
    display = SerialD::ON;
}

void Serial::gen_crc8ccit_table(){
    uint8_t val;
    
    crc8ccit_table = new uint8_t[256];
    for(uint16_t cnt = 0;cnt < 256;cnt++){
        val = (uint8_t)cnt;
        for(uint16_t cnt2 = 0;cnt2 < 8;cnt2++){
            if(val & 0x80){
                val <<= 1;
                val ^= crc8ccit_poly;
            }else{
                val <<= 1;
            }
        }
        crc8ccit_table[cnt] = val;
    }
}

uint8_t Serial::calc_crc8ccit(uint8_t *data,uint8_t len){
    uint8_t crc8ccit = 0x00;
    
    for(uint16_t cnt = 0;cnt < len;cnt++){
        crc8ccit = crc8ccit_table[crc8ccit ^ data[cnt]];
    }
    return crc8ccit;
}


int Serial::read_trush(){
    int len,finishf=1;
    long count = 0;
    int offset=0;
    while(finishf) {
        len = read(fd, buf, sizeof(buf));
        for(int ii = 0; ii < len; ii++) {
            if(display){std::cout << buf[ii] << std::flush;}
            if(buf[ii]=='\n'){
                finishf = 0;
            }   
        }
    }
    if(display){std::cout << std::endl;}
    return 0;
}

int Serial::read_get(uint8_t *buf8t,int len){
    int leng = 0,finishf=1;
    while(leng<len){
        if(leng<1){
            leng = read(fd, buf8t, 1);
        }else{
            leng += read(fd, buf8t+leng, len-leng);
        }
    }
    if(buf8t[len-1] != calc_crc8ccit(buf8t,len-1)){
        if(display){std::cout << "\tfalse" << std::endl;}
        return 0;
    }
    tcflush(fd,TCIFLUSH);
    for(int ii = 0; ii < leng; ii++) {
        if(display){std::cout << std::setfill('0') << std::setw(2) << std::hex << (int)buf8t[ii] << std::flush;}            
    }
    if(display){std::cout << "\ttrue" << std::endl;}
    return 1;
}

int Serial::read_get(uint8_t *buf8t,int len,uint8_t headbyte){
    int leng = 0,finishf=1;
    while(leng<len){
        if(leng<1){
            leng = read(fd, buf8t, 1);
            if(buf8t[0]!=headbyte){
                leng = 0;
            }
        }else{
            leng += read(fd, buf8t+leng, len-leng);
        }
    }
    tcflush(fd,TCIFLUSH);
    if(display){
        for(int ii = 0; ii < leng; ii++) {
            std::cout << std::setfill('0') << std::setw(2) << std::hex << (int)buf8t[ii] << std::flush;            
        }
    }
    return 1;
}

int Serial::read_get(uint8_t *buf8t,int len,int timeoutms){
    int leng = 0,finishf=1;
    int count = 0;
    gettimeofday(&init_time, NULL);
    while(leng<len){  
        if(leng<1){
            leng = read(fd, buf8t, 1);
        }else{
            leng += read(fd, buf8t+leng, len-leng);
        }
        gettimeofday(&end_time, NULL);
        if(((int)((end_time.tv_usec - init_time.tv_usec)/1000) - timeoutms) > 0){return 0;}//タイムアウト
    }
    if(buf8t[len-1] != calc_crc8ccit(buf8t,len-1)){
        if(display){std::cout << "\tfalse" << std::endl;}
        return 0;
    }
    tcflush(fd,TCIFLUSH);
    if(display){
        for(int ii = 0; ii < leng; ii++) {
            std::cout << std::setfill('0') << std::setw(2) << std::hex << (int)buf8t[ii] << std::flush;            
        }
        std::cout << "\ttrue" << std::endl;
    }
    return 1;
}


int Serial::write_string(std::string str){
    ssize_t ret = 0;
    int num = str.size();
    const char* buff = str.c_str();
    std::string strclear = "\r\n";
    const char* buffcl = strclear.c_str();
    ret = write(fd, buffcl, 2);
    ret = write(fd, buff, num);
    ret = write(fd, buffcl, 2);
    return ret;
}

int Serial::write_raw(uint8_t *buf8t,int len){
    int ret = 0;
    ret = write(fd,buf8t,len);
    return ret;
}

int Serial::write_wcrc(uint8_t *buf8t,int len){
    int ret = 0;
    buf8t[len-1] = calc_crc8ccit(buf8t,len-1);
    ret = write(fd,buf8t,len);
    return ret;
}

int Serial::close_s(){
    tcsetattr(fd, TCSANOW, &oldtio);
    close(fd);
    return 0;
}

Serial::~Serial(){
    close_s();
}

#if defined(SERIAL_IS_MAIN)

int main(int argc,char *argv[]){
    std::string devname = "/dev/ttyS11";
    int boardrate = B115200;

    std::string *argvst = new std::string[argc];
    int jj = 0;
    for(int ii=1;ii<argc;ii++){
        argvst[ii] = argv[ii];
        if(argvst[ii].substr(0,4)=="/dev"){
            devname = argvst[ii];
        }
        if(argvst[ii]=="B9600"){
            boardrate = B9600;
        }
        if(argvst[ii]=="B19200"){
            boardrate = B19200;
        }
        if(argvst[ii]=="B38400"){
            boardrate = B38400;
        }
        if(argvst[ii]=="B57600"){
            boardrate = B57600;
        }
        if(argvst[ii]=="B115200"){
            boardrate = B115200;
        }
        if(argvst[ii]=="B230400"){
            boardrate = B230400;
        }
        if(argvst[ii]=="B500000"){
            boardrate = B500000;
        }
        if(argvst[ii]=="B1000000"){
            boardrate = B1000000;
        }
    }
    
    Serial *ser = new Serial(boardrate,devname);
    keyboard ky;

    int motornum = 7;
    int sendbyte = 1+2*motornum+1;
    int readbyte = 3;
    
    
    uint8_t read8t[readbyte];
    uint8_t send8t[sendbyte];
    int readret = 0;
    send8t[0] = 0x0b;
    for(int ii=1;ii<sendbyte-1;ii++){
        send8t[ii] = 0x00;
    }
    struct timeval start_time,end_time;   //時間
    int data[motornum];
    double omega = 3.0;
    double time_second = 0.0;
    ser->setdisplay();
    gettimeofday(&start_time, NULL);
    while(1){
        gettimeofday(&end_time, NULL);
        time_second = (double)(end_time.tv_sec - start_time.tv_sec)+(double)(end_time.tv_usec - start_time.tv_usec)*0.000001;
        data[0] = (int)(120.0*sin(omega*time_second));
        data[1] = (int)(120.0*cos(omega*time_second));
        data[2] = (int)(120.0*sin(0.5*omega*time_second));
        data[3] = (int)(120.0*cos(0.5*omega*time_second));
        data[4] = (int)(120.0*sin(1.5*omega*time_second));
        data[5] = (int)(120.0*cos(1.5*omega*time_second));
        data[6] = (int)(120.0*sin(2.0*omega*time_second));
        for(int ii=0;ii<(sizeof(data)/sizeof(int));ii++){
            send8t[2*ii+1] = ((uint16_t)data[ii] >> 8 ) & 0xFF; 
            send8t[2*ii+2] = (uint16_t)data[ii] & 0xFF;
        }
        ser->write_wcrc(send8t,sizeof(send8t));
        readret = ser->read_get(read8t,sizeof(read8t),5);
        if(0){
            std::cout << std::dec << data << "\t" <<std::setfill('0') << std::setw(2) << std::hex << (int)send8t[1] << (int)send8t[2]<< "\t" <<std::flush;
        }
        //if(readret==0){std::cout << "time out" << std::endl;}
        if(ky.kbhit()){
            break;
        }
    }
    delete ser;
    return 0;
}
#endif