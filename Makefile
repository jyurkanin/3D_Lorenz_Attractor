

CFLAGS= -g -D__LINUX_ASLA -Wall
LIBS= -lm -lpthread -lasound -lX11

all:
	g++ -o ddd main.cpp $(CFLAGS) $(LIBS)
