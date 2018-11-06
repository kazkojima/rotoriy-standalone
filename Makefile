obj= main.o RotorIy.o sincosf.o 
src= $(obj:.cpp=.o)

CPU=$(shell uname -m)
ifeq "$(CPU)" "i686"
ARCH_CPPFLAGS= -msse3 -mfpmath=sse
endif

CPFLAGS= -g -O3 -Wall -ffast-math $(ARCH_CPPFLAGS)
CP= g++

all: rotoriy

rotoriy: $(obj)
	$(CP) $(CPFLAGS) -o $@ $(obj)

%.o: %.cpp
	$(CP) $(CPFLAGS) -c -o $@ $<

RotorIy.o: RotorIy.h
main.o: RotorIy.h

clean:
	rm -f rotoriy *~ *.o
