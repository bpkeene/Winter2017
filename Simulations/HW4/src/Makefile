# Makefile for main.cpp,
# HW4 MENG341 / MolecularDynamics
#
# Compiler: g++

OBJS = main.o simulation.o box.o atom.o \
	   random_mars.o

CC = g++
DEBUG = -g
CFLAGS = -Wall -Wpedantic -c $(DEBUG)
LFLAGS = -Wall -Wpedantic $(DEBUG)

MolecularDynamics : $(OBJS)
	$(CC) $(LFLAGS) $(OBJS) -o MolecularDynamics

main.o: simulation.h 
	$(CC) $(CFLAGS) main.cpp

simulation.o: simulation.h \
	box.h atom.h random_mars.h
	$(CC) $(CFLAGS) simulation.cpp

box.o: box.h atom.h
	$(CC) $(CFLAGS) box.cpp

atom.o: atom.h
	$(CC) $(CFLAGS) atom.cpp

random_mars.o: random_mars.h
	$(CC) $(CFLAGS) random_mars.cpp


clean: 
	rm -f *.o *.dat *.xyz *~ MolecularDynamics
