#include <iostream>
#include <sstream>
#include <fstream>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/ioctl.h>
#include <fcntl.h>
#include <termios.h>
#include <unistd.h>

#include "keyboard.h"

keyboard::keyboard(){tcgetattr(0,&init_tio);}
keyboard::~keyboard(){tcsetattr(0,TCSANOW,&init_tio);}

void keyboard::start_inturkey(){tcgetattr(0,&init_tio);}
void keyboard::finish_inturkey(){tcsetattr(0,TCSANOW,&init_tio);}

int keyboard::kbhit(){
    struct termios tio;
    struct timeval tv;
    fd_set rfds;
    memcpy(&tio,&init_tio,sizeof(struct termios));
    tio.c_lflag &= ~(ICANON);
    tcsetattr(0,TCSANOW,&tio);
    FD_ZERO(&rfds);
    FD_SET(0,&rfds);
    tv.tv_usec = 0;
    tv.tv_sec  = 0;
    select(1,&rfds,NULL,NULL,&tv);
    tcsetattr(0,TCSANOW,&init_tio);
    return (FD_ISSET(0,&rfds)?1:0);
}
    
int keyboard::getkey(){
    if(kbhit()){
        int ch;
        struct termios tio;
        memcpy(&tio,&init_tio,sizeof(struct termios));
        tio.c_lflag &= ~(ICANON | ECHO);
        tcsetattr(STDIN_FILENO,TCSANOW,&tio);
        read(0,&ch,1);
        tcsetattr(0,TCSANOW,&init_tio);
        return ch;
    }
    return 0;
}

#if defined(KEY_IS_MAIN)
int main(){
    keyboard *key;
    key = new keyboard; 
    char ch;
    while(!key->kbhit()){
      std::cout << "input key------"<<std::endl;  
    }
    key->finish_inturkey();
    while(1){
        std::cout << "hogehoge"<<std::endl; 
        ch = getchar();
        std::cout << "input key is " << ch << std::endl;
    }
    return 0;
}
#endif