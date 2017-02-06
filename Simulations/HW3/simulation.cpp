#include "simulation.h"
#include <iostream>
#include <math.h>
#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>

// declare our constructor
Simulation::Simulation(int _numberOfAtoms, double _density, 
                   double _sigma, double _epsilon, int _seed) :
            numberOfAtoms(_numberOfAtoms), density(_density), 
            sigma(_sigma), epsilon(_epsilon)
{
// and stuff to do here
    prng = new RanMars(_seed);
    box = new Box(numberOfAtoms, density, sigma);
    step = 0;
    distances = std::vector< std::vector<double> > (numberOfAtoms, std::vector<double> (numberOfAtoms));
    ijDistance = std::vector<double> (numberOfAtoms, 0.0);
};


void Simulation::initializeAtoms() {
    // initialize the atoms on the lattice in the box;
    // what should be the inter-atomic distance?

    // initialize the vector of atoms
    atoms = std::vector<Atom> (numberOfAtoms);

    // and tell box to place the atoms in the box & compute distances subject to the geometry
    box->initializeAtoms(atoms,distances);
    std::cout << "returned from box->initializeAtoms" << std::endl;
}; 

// our equilibration steps
void Simulation::run(int nsteps) {
    // for a given number of steps..
    accepted = 0.0;
    total = 0.0;

    // and we'll do block averages as well for sampling for good alpha
    double acceptedInterval = 0.0;
    double totalInterval = 0.0;

    alpha = 1.0; // initialize to some value, we'll say 1.0
    double r1, r2, r3; // our three random numbers
    int idx; // our randomly selected index - the atom which we will move
    
    double pre_move_pot = 0.0;
    double post_move_pot = 0.0;
  
    // every 5% of moves, change alpha by some factor if local acceptance is not 40 < x < 60 %
    int changeAlphaEvery = (int) floor(0.05 * nsteps);

    int nAtomsIdx = atoms.size() - 1;
    // get the current potential energy of the system
    ComputeTotalEnergy();

    for (int i = 0; i < nsteps; i++) {

        // randomly select an atom to move
        idx = lrint(prng->uniform() * nAtomsIdx);
        
        // save its coordinates
        atom[idx].saveCoordinates();
        
        // compute the energy of its current coordinates
        pre_move_pot = ComputeAtomPotential(idx);
        
        // translate the atom, applying PBC
        box->Translate(atom[idx], alpha, prng); 

        // compute delta E
        post_move_pot = ComputeAtomPotential(idx);

        // accept or reject the move
        acceptMove = acceptOrReject(pre_move_pot, post_move_pot);

        if (acceptMove) {
            acceptedInterval += 1.0;
            accepted += 1.0;
            // update the distances in our distance matrix
            updateDistances(idx);

            totalPE += (post_move_pot - pre_move_pot);

        } else {
            // set the coordinates of this atom back to its old coordinates
            atoms[idx].resetCoordinates();
        };

        // update the total move counters, for this interval and for the simulation
        totalInterval += 1.0;
        total += 1.0;

        if ( ( (i % changeAlphaEvery) == 0) and (i != 0)) {
            double successRatioInterval = acceptedInterval / totalInterval;
            if (successRatioInterval >= 0.55) {
                // we can increase alpha a bit to take larger steps
                alpha *= 1.05;
            } else if (successRatioInterval <= 0.45) {
                // we should take smaller steps..
                alpha *= 0.95;
            };
            std::cout << "At step " << i << " of equilibration " << "with acceptance ratio " << successRatioInterval << std::endl;
        };
    };
};


// our production steps
void Simulation::run(int nsteps,int printEvery) {
    


};

bool Simulation::acceptOrReject(double deltaE) {

    return true;

};


void Simulation::ComputeTotalEnergy() {
    // zero our totalPE variable
    totalPE = 0.0;
    for (int i = 0; i < (distances[0].size() - 1); i++) {
        for (int j = i+1; j < distances[0].size(); j++) {
            totalPE += LJPotential(i,j);
        };
    };
};

double Simulation::LJPotential(Atom &atom1, Atom &atom2) {
    double rij = distances[i][j];
    double sigma = atoms[i].getSigma();
    double epsilon = atoms[i].getEpsilon();

    double ijEnergy = 0.0;
    ijEnergy = 4.0 * epsilon ; 

    return ljEnergy;

};

void Simulation::printConfig(std::string name, int step) {
    std::ostringstream stringStream;
    stringStream.flush();
    stringStream.str("");
    stringStream << name << "_step_" << step << ".xyz";

    std::string xyzFileName = stringStream.str();

    std::ofstream xyzFile(xyzFileName.c_str(), std::ios::out);
    std::vector<double> boxDim = box->getDim();
    xyzFile << "# " << atoms.size() << " Monte Carlo simulation\n" << std::endl;
    xyzFile << "#\n#\n" << "# format: xpos ypos zpos" << std::endl;
    
    xyzFile << "# box dimensions:" << std::endl;
    xyzFile << "# " << boxDim[0] << " " << boxDim[1] << " " << boxDim[2] << std::endl;
    
    xyzFile << atoms.size() << "\n" <<  std::endl;

    std::vector<double> theseCoords = std::vector<double> (3,0.0);
    for (unsigned int i = 0; i < atoms.size(); i++) {
        theseCoords = atoms[i].getCoordinates();
        xyzFile << i << " " << theseCoords[0] << " " << theseCoords[1] << " " << theseCoords[2] << std::endl;
    };


};





