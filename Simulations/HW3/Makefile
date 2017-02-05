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

main.o: simulation.h 
	$(CC) $(CFLAGS) main.cpp

simulation.o: simulation.h \
	atom.h box.h random_mars.h
	$(CC) $(CFLAGS) simulation.cpp

atom.o: atom.h
	$(CC) $(CFLAGS) atom.cpp

box.o: box.h atom.h
	$(CC) $(CFLAGS) box.cpp

random_mars.o: random_mars.h
	$(CC) $(CFLAGS) random_mars.cpp


# these are correct; the ones above are not (necessarily)
#simulator.o: simulator.cpp random_mars.h
#	$(CC) $(CFLAGS) simulator.cpp

#random_mars.o : random_mars.h random_mars.cpp
#	$(CC) $(CFLAGS) random_mars.cpp

clean: 
	rm -f *.o *~ MonteCarlo






