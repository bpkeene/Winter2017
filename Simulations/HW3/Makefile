# Makefile for main.cpp
#
# Compiler: g++

OBJS = main.o simulation.o atom.o \
	   random_mars.o

CC = g++
DEBUG = -g
CFLAGS = -Wall -c $(DEBUG)
LFLAGS = -Wall $(DEBUG)

MonteCarlo : $(OBJS)
	$(CC) $(LFLAGS) $(OBJS) -o MonteCarlo

main.o: main.cpp simulation.h 
	$(CC) $(CFLAGS) main.cpp

simulation.o: simulation.h simulation.cpp


# these are correct; the ones above are not (necessarily)
#simulator.o: simulator.cpp random_mars.h
#	$(CC) $(CFLAGS) simulator.cpp

#random_mars.o : random_mars.h random_mars.cpp
#	$(CC) $(CFLAGS) random_mars.cpp

clean: 
	rm -f *.o *~ MonteCarlo






