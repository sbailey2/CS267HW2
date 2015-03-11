#
# Hopper - NERSC 
#
# Portland Group Compilers PGI are loaded by default; for other compilers please check the module list
#
CC = CC
MPCC = CC
OPENMP = -mp
CFLAGS = -O3
LIBS =


#TARGETS = serial autograder pthreads mpi
TARGETS = autograder mpi

all:	$(TARGETS)

serial: serial.o common.o
	$(CC) -o $@ $(LIBS) serial.o common.o
autograder: autograder.o common.o
	$(MPCC) -o $@ $(LIBS) autograder.o common.o
#pthreads: pthreads.o pthread_barrier.o common.o
#	$(CC) -o $@ $(LIBS) -lpthread pthreads.o pthread_barrier.o common.o
pthreads: pthreads.o common.o
	$(CC) -o $@ $(LIBS) -lpthread pthreads.o common.o 
openmp: openmp.o common.o
	$(CC) -o $@ $(LIBS) $(OPENMP) openmp.o common.o
mpi: mpi.o common2.o
	$(MPCC) -o $@ $(LIBS) $(MPILIBS) mpi.o common2.o

autograder.o: autograder.cpp common.h
	$(MPCC) -c $(CFLAGS) autograder.cpp
openmp.o: openmp.cpp common.h
	$(CC) -c $(OPENMP) $(CFLAGS) openmp.cpp
serial.o: serial.cpp common.h
	$(CC) -c $(CFLAGS) serial.cpp
pthreads.o: pthreads.cpp common.h
	$(CC) -c $(CFLAGS) pthreads.cpp
mpi.o: mpi.cpp common.h
	$(MPCC) -c $(CFLAGS) mpi.cpp
common.o: common.cpp common.h
	$(MPCC) -c $(CFLAGS) common.cpp
common2.o: common2.cpp common2.h
	$(MPCC) -c $(CFLAGS) common2.cpp
#pthread_barrier.o: pthread_barrier.c pthread_barrier.h 
#	$(CC) -c pthread_barrier.c

clean:
	rm -f *.o $(TARGETS) *.stdout *.txt
