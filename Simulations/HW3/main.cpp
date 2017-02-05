#include "simulation.h" // probably to run a simulation
#include <iostream> // for the std::cout
#include <sstream> // for std::stringstream
#include <string>
#include <cstdio>
#include <fstream>
#include <stdlib.h> 
#include <stdio.h>
#include <string.h>

int main(int argc, char *argv[]) {

    std::string filename;

    if (argc != 2) {
        std::cout << "No input file provided, or too many arguments.  Exiting." << std::endl;
        return(1);
    }; 
    
    if ((argv[1] == NULL)) {
        return(1);
    };

    std::ifstream fileStream(argv[1]);

    FILE *input;
    if ( NULL == (input = fopen(argv[1],"r"))) {
        fprintf(stdout, "Error reading input file. Could not find.\n");
            return(1);
    };

    input = fopen(argv[1],"r");

    /* 
     * We look for the following parameters in the input file:
     * - a density (reduced)
     * - # of atoms
     * - LJ sigma parameter
     * - LJ epsilon parameter
     *
     */

    double density; // reduced units
    int nAtoms; // integer
    double sigma; // regular units
    double epsilon; // regular units
    int seed; // seed for Marsaglia PRNG
    int equilibration; // number of steps
    int production; // number of steps
    int printXYZ; // every number of steps
    char simulationNameC[256]; // fscanf reads chars
    std::string simulationName; // recast as string
    // buffer by which to move forward;
    char tt[80];
    
    // assume the first line contains the density
    // this code for reading input was taken from ljmd.cpp
    fscanf(input, "%lf", &density); fgets(tt, 80, input);
    fscanf(input, "%d", &nAtoms); fgets(tt, 80, input);
    fscanf(input, "%lf", &sigma); fgets(tt, 80, input);
    fscanf(input, "%lf", &epsilon); fgets(tt, 80, input);
    fscanf(input, "%d", &seed); fgets(tt,80,input);
    fscanf(input, "%d", &equilibration); fgets(tt, 80, input);
    fscanf(input, "%d", &production); fgets(tt, 80, input);
    fscanf(input, "%d", &printXYZ); fgets(tt, 80, input);
    fscanf(input, "%s", simulationNameC); fgets(tt, 80, input);
    fclose(input);
    // end of code taken from ljmd.cpp

    // we prefer working with strings, so do a type change
    simulationName = "";
    simulationName += simulationNameC;

    
    std::cout << "value of density: " << density << std::endl;
    std::cout << "value of nAtoms: " << nAtoms << std::endl;
    std::cout << "value of sigma: " << sigma << std::endl;
    std::cout << "value of epsilon: " << epsilon << std::endl;
    std::cout << "value of seed: " << seed << std::endl;
    std::cout << "value of equilibration: " << equilibration << std::endl;
    std::cout << "value of production: " << production << std::endl;
    std::cout << "value of printXYZ: " << printXYZ << std::endl;
    //std::cout << "value of simulation name: " << simulationName << std::endl;
    

    
    // define the simulation state
    Simulation simulation(nAtoms, density, sigma, epsilon, seed);

    // initialize the atoms on a lattice
    simulation.initializeAtoms();
    
    // print the initial configuration; simply pass in the simname and step count (0)
    simulation.printConfig(simulationName, 0);
    
    /*
    // run the equilibration steps
    simulation.run(equilibration);

    // run the production steps
    simulation.run(production,printXYZ,trajectory);
    */
    return 0;
};













