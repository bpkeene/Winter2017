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
    box = new Box(numberOfAtoms, sigma, density);
    step = 0;
    distances = std::vector< std::vector<double> > (numberOfAtoms, std::vector<double> (numberOfAtoms));
};


void Simulation::initializeAtoms() {
    // initialize the atoms on the lattice in the box;
    // what should be the inter-atomic distance?
    // --- if hard sphere potential, initialize s.t.
    double increment;
    double initial_x_pos, initial_y_pos, initial_z_pos;
    double xpos, ypos, zpos;
    atoms = std::vector<Atom> (numberOfAtoms);
    
    // cast number of atoms as double
    double nAtoms = floor(numberOfAtoms);

    // atoms per length is N^(1/3); put as ceiling, to make 
    // sure we can fit all the atoms
    double atomsPerLengthD = ceil(pow(nAtoms, (1.0/3.0))); 

    // get the dimensions of our box
    std::vector<double> maxDimensions = box->getDim();


    // divide x-dim by atoms per length to get the increment
    // since it is cubic, this increment will be the same in all directions
    increment = maxDimensions[0] / atomsPerLengthD;
    int atomsPerLength = lrint(atomsPerLengthD);

    // we evenly distribute the atoms on a cubic lattice;
    initial_x_pos = 0.0;
    initial_y_pos = 0.0;
    initial_z_pos = 0.0;

    xpos = ypos = zpos = 0.0;

    int atomIndex = 0;

    while (atomIndex < nAtoms) {
        for (int i = 0; i < atomsPerLength; i++) {
            zpos = initial_z_pos + i * increment;
                for (int j = 0; j < atomsPerLength; j++) {
                    ypos = initial_y_pos + j * increment;
                    xpos = initial_x_pos;
                    for (int k = 0; k < atomsPerLength; k++) {
                        atoms[atomIndex].setCoordinates(xpos,ypos,zpos);
                        xpos += increment;
                        atomIndex += 1;
                    };
                };
        };
    };

    // we also may as well initialize the matrix containing the distances between all atoms
    for (unsigned int m = 0; m < distances[0].size(); m++) {
        for (unsigned int n = 0; n < distances[0].size(); n++) {
            if (n <= m) {
                distances[m][n] = 0.0;
            } else {
                distances[m][n] = box->computeDistance(atoms[m],atoms[n]);
            };
        };
    };
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





