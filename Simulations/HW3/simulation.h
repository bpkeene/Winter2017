#ifndef SIMULATION_H
#define SIMULATION_H

#include "atom.h"
#include "box.h"
#include "random_mars.h"
#include <vector>
#include <math.h>
#include <string>

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
        int step; // current step in this run; initialized to 0, and set to 0 on calling 'run()'
        int printEvery; // print every number of steps

        double accepted, total; // accepted and total number of steps
        bool acceptMove; // accept the move or not
        double alpha; // factor for the random numbers, determining the size of our translation

        double totalPE; // our total potential energy within the system

        // we'll have a box of type box
        Box *box;

        // a vector of type Atom; this is where we have all our atoms
        std::vector<Atom> atoms;

        // and we'll have a vector of all i-j distances; this will be updated
        // on moving
        std::vector< std::vector<double> > distances;
     
        // some temporary vector of length N containing some distances between atom i and other
        // atoms j
        std::vector<double> ijDistance;

        // method for determining whether to accept or a reject a move
        bool acceptOrReject(double deltaE);

        RanMars *prng;

        void Translate(int nSteps);

        void ComputeTotalEnergy();
        
        double LJPotential(Atom &, Atom &);

        void saveDistances(Atom &);

    public:

        Simulation(int _numberOfAtoms, double _density, 
                   double _sigma, double _epsilon, int _seed); 

        // initialize the atoms on a lattice
        void initializeAtoms();
        
        // runs for number of steps
        void run(int);

        // runs for number of steps, and prints config every number of steps
        void run(int,int);

        // prints the configuration
        void printConfig(std::string name, int step);
};


#endif // end of simulator.h
