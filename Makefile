#
# Hopper - NERSC 
#
# Portland Group Compilers PGI are loaded by default; for other compilers please check the module list
#
CC = CC
MPCC = CC
OPENMP = -mp
CFLAGS = -g -std=c++11
LIBS = -lstdc++

TARGETS = serial autograder
#TARGETS = serial pthreads openmp mpi autograder

all:	$(TARGETS)

serial: serial.o common.o
	$(CC) -o $@ serial.o common.o $(LIBS)
autograder: autograder.o common.o 
	$(CC) -o $@ autograder.o common.o $(LIBS)
#pthreads: pthreads.o common.o
#	$(CC) -o $@ $(LIBS) -lpthread pthreads.o common.o
#openmp: openmp.o common.o
#	$(CC) -o $@ $(LIBS) $(OPENMP) openmp.o common.o
#mpi: mpi.o common.o
#	$(MPCC) -o $@ $(LIBS) $(MPILIBS) mpi.o common.o

autograder.o: autograder.cpp common.h
	$(CC) -c $(CFLAGS) autograder.cpp
#openmp.o: openmp.cpp common.h
#	$(CC) -c $(OPENMP) $(CFLAGS) openmp.cpp
serial.o: serial.cpp common.h
	$(CC) -c $(CFLAGS) serial.cpp
#pthreads.o: pthreads.cpp common.h
#	$(CC) -c $(CFLAGS) pthreads.cpp
#mpi.o: mpi.cpp common.h
#	$(MPCC) -c $(CFLAGS) mpi.cpp
common.o: common.cpp common.h
	$(CC) -c $(CFLAGS) common.cpp

clean:
	rm -f *.o $(TARGETS) *.stdout *.txt
