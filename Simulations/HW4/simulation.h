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

// class Simulation
class Simulation {
    
    private:
        int numberOfAtoms;

        double mass;
        double dt; // our integration timestep
        double dt2m; // dt/2m, used in integrating our velocities
        double dt2; // dt/2, used in integrating our positions
        double density;
        double sigma; // put these in the atoms class later
        double sig3;
        double sig6;
        double sig12;

        double epsilon; // put this in the atoms class later
        double Tstar;

        double tail; // our tail correction (constant, once rcut is set)
        
        int seed; // seed for PRNG; we use this to initialize the velocities
        int step; // current step in this run; initialized to 0, and set to 0 on calling 'run()'
        int printEvery; // print every number of steps
        double rcut; // nominally set to 3 * sigma;
        double rcut3; // rcut ^ 3
        double rcut9; // rcut ^ 9
        double total; // total number of moves

        std::string name;
       
        std::vector<double> deltaR2; // for calculating the diffusionCoeff;

        double totalPE; // our total potential energy within the system
        double totalKE; // total kinetic energy (3/2 * N k T)
        double totalKE_n; // same quantity, but 3/2 k T

        // we'll have a box of type box
        Box *box;

        // a vector of type Atom; this is where we have all our atoms
        std::vector<Atom> atoms;

        RanMars *prng;

        void ComputeTotalPotentialEnergy();
        
        double LJPotential(int, int);
        
        void LJForce(int, int);

        void ComputeTotalForce();

        double ComputeAtomPotential(int);
   
        void velocityHalfStep();

        void positionStep();

        void zeroForces();

        void calculateForces();

        void calculateAvgDeltaR2();

    public:

        // our constructor
        Simulation(int _numberOfAtoms, double _density, double _dt, double _mass,
                   double _sigma, double _epsilon, double _Tstar, int _seed,
                   std::string _name); 

        // initialize the atoms on a lattice
        void initializeAtomsPositions();
       
        // initialize the atoms velocities, sampling from the boltzmann distribution
        void initializeAtomsVelocities();

        // runs for number of steps
        void run(int, int, bool);

        // prints the configuration
        void printConfig(std::string name, int step);
};


#endif /* SIMULATION_H */
