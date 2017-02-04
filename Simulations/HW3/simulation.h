#ifndef SIMULATION_H
#define SIMULATION_H






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
        int step; // current step in this run
        int printEvery; // print every number of steps
    
    public:
        Simulation(int _numberOfAtoms, double _density, 
                   double _sigma, double _epsilon) :
            numberOfAtoms(_numberOfAtoms), density(_density), 
            sigma(_sigma), epsilon(_epsilon) {
            };

        // runs for number of steps
        void run(int);

        // runs for number of steps, and prints config every number of steps
        void run(int,int);

        // prints the configuration
        void printConfig();
};


#endif // end of simulator.h
