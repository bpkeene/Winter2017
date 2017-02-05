#include "simulation.h"


// declare our constructor
Simulation::Simulation(int _numberOfAtoms, double _density, 
                   double _sigma, double _epsilon, int _seed) :
            numberOfAtoms(_numberOfAtoms), density(_density), 
            sigma(_sigma), epsilon(_epsilon)
{
// and stuff to do here
    prng = new RanMars(_seed);
    box = Box(_numberOfAtoms, _sigma, _density);

};


void Simulation::run(int nsteps) {




};

void Simulation::run(int nsteps,int printEvery) {



};

void Simulation::printConfig() {




};





