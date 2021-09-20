#ifndef KBHIT_H
#define KBHIT_H

class keyboard{
    protected:
        struct termios init_tio;
    public:
        keyboard();
        ~keyboard();
        void start_inturkey();
        void finish_inturkey();
        int kbhit();
        int getkey();
};

#endif