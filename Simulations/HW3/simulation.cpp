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



};



void Simulation::run(int nsteps,int printEvery) {



};

bool Simulation::acceptOrReject(double deltaE) {

    return true;

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





