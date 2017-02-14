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

    double density; // reduced units
    int nAtoms; // integer
    double mass; // amu? idk
    double sigma; // regular units
    double epsilon; // regular units
    double Tstar; // dimensionless temperature
    int seed; // seed for Marsaglia PRNG
    double timestep; // timestep to be taken; coefficient to LJ time
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
    fscanf(input, "%lf", &mass); fgets(tt, 80, input);
    fscanf(input, "%lf", &sigma); fgets(tt, 80, input);
    fscanf(input, "%lf", &epsilon); fgets(tt, 80, input);
    fscanf(input, "%lf", &Tstar); fgets(tt,80,input);
    fscanf(input, "%d", &seed); fgets(tt,80,input);
    fscanf(input, "%lf", &timestep); fgets(tt,80,input);
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
    std::cout << "value of mass: " << mass << std::endl;
    std::cout << "value of sigma: " << sigma << std::endl;
    std::cout << "value of epsilon: " << epsilon << std::endl;
    std::cout << "value of Tstar: " << Tstar << std::endl;
    std::cout << "value of seed: " << seed << std::endl;
    std::cout << "value of timestep: " << timestep << std::endl;
    std::cout << "value of equilibration: " << equilibration << std::endl;
    std::cout << "value of production: " << production << std::endl;
    std::cout << "value of printXYZ: " << printXYZ << std::endl;
    std::cout << "value of simulation name: " << simulationName << std::endl;
    
    // define the simulation state
    Simulation simulation(nAtoms, density, timestep, mass, sigma, epsilon, Tstar, seed, simulationName);

    // initialize the atoms on a lattice
    simulation.initializeAtomsPositions();
 
    //initialize the atoms with some velocities from the Boltzmann distribution;
    //note that these will then also be scaled s.t. momentum x, y, z is zero
    simulation.initializeAtomsVelocities();
    
    // we are doing molecular dynamics; give them velocities as well,
    // sampled from the Boltzmann distribution
    // print the initial configuration; simply pass in the simname and step count (0)
    
    simulation.printConfig(simulationName, 0);
    
    bool isProduction = false;

    //TODO alert: melding equilibration and production
    equilibration += production;
    // run the equilibration steps
    simulation.run(equilibration,printXYZ,isProduction);
    
    isProduction = true;
    std::cout << "completed equilibration!" << std::endl;
    
    // run the production steps
    //simulation.run(production,printXYZ,isProduction);
    
    simulation.printConfig(simulationName,production);

    return 0;
}













