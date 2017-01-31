# Makefile for simulator.cpp with 
# dependencies on random_mars.h
#
# Compiler: g++
#

OBJS = simulator.o random_mars.o
CC = g++
DEBUG = -g
CFLAGS = -Wall -c $(DEBUG)
LFLAGS = -Wall $(DEBUG)

MonteCarlo : $(OBJS)
	$(CC) $(LFLAGS) $(OBJS) -lfftw3 -lm -o MonteCarlo

simulator.o: simulator.cpp random_mars.h
	$(CC) $(CFLAGS) simulator.cpp

random_mars.o : random_mars.h random_mars.cpp
	$(CC) $(CFLAGS) random_mars.cpp

clean: 
	rm -f *.o *.xyz *.rdf *.fes *.sf *~ MonteCarlo






