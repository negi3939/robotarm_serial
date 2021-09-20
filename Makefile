ifdef target
	TARGET=$(target)
else
	TARGET=serial
endif

ifeq ($(TARGET),key)
	SOURCE_MAIN = keyboard.cpp
	CXXFLAGS = -DKEY_IS_MAIN
endif

ifeq ($(TARGET),serial)
	SOURCE_MAIN = serial.cpp
	SOURCE_SUB = keyboard.cpp
	CXXFLAGS = -DSERIAL_IS_MAIN
endif


PROGRAM = $(SOURCE_MAIN:%.cpp=%.exe)
SUBOBJ = $(SOURCE_SUB:%.cpp=%.o)

DIRX = /usr/X11R6/lib

CXXFLAGS += -fpermissive
LDFLAGS	 = -L "$(DIRX)" -lm
CXX = g++-7

all: $(PROGRAM)

%.out: %.o $(SUBOBJ)
	$(CXX) -o $@ $^ $(LDFLAGS) -w
%.exe: %.o $(SUBOBJ)
	$(CXX) -o $@ $^ $(LDFLAGS) -w
%.o : %.cpp
	$(CXX) -o $@ $< -c $(CXXFLAGS) -w
clean:
	rm -f *.o $(PROGRAM)