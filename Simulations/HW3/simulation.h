#ifndef SIMULATION_H
#define SIMULATION_H

#include "atom.h"
#include "box.h"
#include "random_mars.h"
#include <vector>
#include <math.h>

/* Simulation.h
 * Holds the main that define our simulation, and runs 
 * it; essentially our 'main', less the input file.
 */

class Simulation {
    
    private:
        int numberOfAtoms;
        double density;
        double sigma; // put these in the atoms class later
        double epsilon; // put this in the atoms class later
        int seed; // seed for PRNG
        int step; // current step in this run
        int printEvery; // print every number of steps

        // we'll have a box of type box
        Box box;

        // a vector of type Atom; this is where we have all our atoms
        std::vector<Atom> atoms;

        // and we'll have a vector of all i-j distances; this will be updated
        // on moving
        std::vector< std::vector<double> > distances;
        
        RanMars *prng;

        void Translate(int nSteps);

    public:
        Simulation(int _numberOfAtoms, double _density, 
                   double _sigma, double _epsilon, int _seed); 

        // runs for number of steps
        void run(int);

        // runs for number of steps, and prints config every number of steps
        void run(int,int);

        // prints the configuration
        void printConfig();
};


#endif // end of simulator.h
