EXEC=sr
CC=gcc
CXX=g++
LD=g++
DBG_CFLAGS= -Wall -D_DEBUG -g
CFLAGS= $(DBG_CFLAGS)
CXXFLAGS=$(DBG_CFLAGS)
LDFLAGS=-L/usr/local/lib -fopenmp
LIBS= -lopencv_core -lopencv_highgui -lopencv_imgproc 
OBJS=img_utils.o backproject.o Demo_SR.o ScSR.o alloc_util.o clust_invert.o l1qp_feasign.o

all: $(OBJS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $^ -o $(EXEC) $(LIBS) 

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -c -o $@ $<

clean:
	rm *.o $(EXEC)
