/* a series of include statements */
#include <math.h>
#include "random_mars.h"
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include "time.h"
#include "stdlib.h"
#include <vector>
#include <fenv.h> 
//#include <fftw3.h>


// global variable declarations
double d0; // the excluded zone (a diameter)
double alpha; // function of d0
int steps;
//double nuValues[8] = {2.0,4.0,5.0,5.5 , 6 , 6.25,6.5,7.0};  
//double KValues[8] = { 1.1,1.3,1.5,1.55,1.6, 1.7 ,1.8,2.0}
// for now, just look at nu = 2
double nuValues[1] = {5.0};
double KValues[1] = {6.0};
double successRatio;
double accepted;
double total;

// for naming files for output
std::ostringstream stringStream;


// we make a class 'Atom' to describe a given hard-sphere atom
// really just need the diameter and coordinates of a given atom
class Atom {
    public:
        double diameter;
        double diameterSq;
        std::vector<double> coords;

        // our constructor
        Atom() {
            diameter = -1.0;
            coords = std::vector<double> (2, 0.0); // initialize array of size 2 to hold coordinates (0.0, 0.0);
            diameterSq = 1.0;
        };

        // function by which to set 'coords'
        void setCoordinates(double x_, double y_) {
            coords[0] = x_;
            coords[1] = y_;
        };

        void setDiameter (double diameter_) {
            diameter = diameter_;
            diameterSq = diameter_ * diameter_;
        };


};

// our simdata class; holds the box dimensions and implements pbc
class Box {
    public:
        double x_dim;
        double y_dim;
        double rn;
        double xpos;
        double ypos;

        double dx, dy, d1,d2;

        Box(double x_dim_, double y_dim_) {
            x_dim = x_dim_;
            y_dim = y_dim_;
        };

        // takes two atoms and computes the distance between their minimum images
        double computeDistance(Atom &atom1, Atom &atom2) {
            dx = atom1.coords[0] - atom2.coords[0];
            dy = atom1.coords[1] - atom2.coords[1];
            
            d1 = fabs(fabs(dx) - x_dim);
            d2 = fabs(dx);
            dx = std::min(d1,d2);

            d1 = fabs(fabs(dy) - y_dim);
            d2 = fabs(dy);
            dy = std::min(d1,d2);
            
            // return the distance and compare to diameterSq;
            return sqrt(dx*dx + dy*dy);
        };    
        
        // modifies the atoms position by some random amount, subject to PBC;
        // note: the old position should be saved /before/ this is called
        void Translate(Atom &atom, double alpha_, RanMars *prng) {
            xpos = atom.coords[0];
            ypos = atom.coords[1];
            
            // our random number
            //double rn = (prng->uniform() * 2.0 ) - 1.0;
            rn = (prng->uniform() * 2.0 ) - 1.0;
            xpos += (alpha_ * rn);

            rn = (prng->uniform() * 2.0 ) - 1.0;
            ypos += (alpha_ * rn);

            // PBC constraints; make sure its still in the square
            // if we move beyond {x || y} == 1, subtract {xdim || ydim}
            // if we move below  {x || y} == 0, add {xdim || ydim}
            if (xpos < 0.0) {
                xpos += x_dim; 
            };
            if (xpos > x_dim) {
                xpos -= x_dim;
            };

            if (ypos < 0.0) {
                ypos += y_dim;
            };

            if (ypos > y_dim) {
                ypos -= y_dim;
            };

            // write our new values to atom.coords
            atom.coords[0] = xpos;
            atom.coords[1] = ypos;
            

        }
};

// a class holding our RDF data; 
class RDF {
    public:
        
        // this object will hold the number of bins M,
        // the width of the bins, specified by our K value
        // and our histogram object, gr
        // bins
        int M;
        double dr;
        double lowerBound;
        double upperBound;
        double d0;

        // our cumulative histogram
        std::vector<double> gr;
        double Kval;

        // a picture of our RDF at one timestep
        std::vector<double> snapshot;

        // our single atom gr modifications of snapshot
        std::vector<double> pre_move_ij;
        std::vector<double> post_move_ij;

        RDF(int M_, double lowerBound_, double upperBound_, double K_) {
            M = M_;
            lowerBound = lowerBound_;
            upperBound = upperBound_;

            // M_ number of bins initialized to zero
            gr = std::vector<double> (M_,0.0);
            // additionally; we can now define our dr
         
            // initialize our snapshot;
            snapshot = std::vector<double> (M_, 0.0);
            // also, initialize our temporary histograms that will be used to update snapshot
            pre_move_ij = std::vector<double> (M_, 0.0);
            post_move_ij = std::vector<double> (M_, 0.0);

            // our Kval is given by ... the K_ that we passed in to the function
            Kval = K_;
            
            d0 = upperBound / Kval;

            // likewise, if our bins start at 0,
            dr = (upperBound_ - lowerBound_) / (double (M_));

        };

        // populate the initial snapshot; note that we do not need to modify distanceMatrix
        void initialStep(const std::vector<std::vector<double> > &distanceMatrix) {
            double dist;
            int binIdx;
            for (unsigned int i = 0; i < (distanceMatrix[0].size() - 1); i++) {
                for (unsigned int j = i+1; j < distanceMatrix[0].size(); j++) {
                    dist = distanceMatrix[i][j];
                    if (dist <= upperBound && dist >= lowerBound) {
                        binIdx = lrint(floor( (dist - lowerBound) / dr));
                        snapshot[binIdx] += 1;
                    }
                }
            }
        };

        void updateSnapshot(std::vector<double> &pre_attempt, std::vector<double> &post_attempt, int idx) {
            // so, we calculate the vector of pre_move ij, calculate post_move_ij, and send the difference to snapshot
            
            // first, send the values contained in pre_move_ij and post_move_ij to zero
            pre_move_ij = std::vector<double> (M, 0.0);
            post_move_ij = std::vector<double> (M, 0.0);
            double dist;
            int binIdx;

            for (unsigned int i = 0; i<pre_attempt.size(); i++) {
                // populate pre_move_ij rdf
                dist = pre_attempt[i];
                if (dist <= upperBound && dist >= lowerBound) {
                    binIdx = lrint(floor( (dist - lowerBound) / dr));
                    pre_move_ij[binIdx] += 1;
                };

                // populate post_move_ij rdf
                dist = post_attempt[i];
                if (dist <= upperBound && dist >= lowerBound) {
                    binIdx = lrint(floor( (dist - lowerBound) / dr));
                    post_move_ij[binIdx] += 1;
                }
            };

            for (unsigned int j = 0; j < snapshot.size(); j++) {
                snapshot[j] += post_move_ij[j] - pre_move_ij[j];
            };
        };

        // add the snapshot of the current system to our cumulative gr histogram
        // ---> gr is simply a sum of the snapshots;
        void calculateRDF() {
            for (unsigned int i = 0; i < gr.size(); i++) {
                gr[i] += snapshot[i];
            };
        }

        // we must normalize the RDF at the end of the production run;
        // just operate on the data stored within the class; easier this way
        void normalizeRDF(int nSteps) {
            double data = 0.0;
            double thisRadius = 0.0;
            double normalization;

            for (unsigned int i = 0; i < gr.size(); i++) {
                data = gr[i];

                // the radius at which we normalize is the radius denoting the /center/ of the histogram bin
                // for i = 0 (the first bin), this is the lowerBound + 0.5 * dr;
                // for i = 1 (the second bin), this is the lowerBound + 0.5*dr + 1 * dr
                // for i = 2 ..................................................+ 2 * dr
                thisRadius = lowerBound + (i * dr) + (0.5 * dr);

                // our normalization is (nAtoms * nAtoms) * $\pi$ * r * dr * (number of steps);
                normalization = (224 * 224) * M_PI * thisRadius * dr * nSteps;

                gr[i] = data / normalization;


            };
        };

        // since the rdf is already here, may as well calculate the corresponding free energy surface;
        // and, when we write the rdf to file, we will write the fes data to file as well
        void calculateFES() {
            // initialize some class members to put in here so that the data doesn't die with the function
        
        };
};

// function declarations 

// Initialize our array of atoms, specifying their coordinates;
// we take our atomlist in, our box, and the distance matrix (to populate it)
void InitializeCoordinates(std::vector<Atom> &, Box &, std::vector<std::vector<double> > &);

// Write a .xyz file of the current configuration;
// take the list of atoms, our value of nu, and the current step 
void PrintConfigurationData(std::vector<Atom> &, double, int);

// takes a vector of atoms, the box, pointer to our PRNG, number of steps, 
// and whether or not we are doing a production run
void MonteCarloTranslate(std::vector<Atom> &, Box &, std::vector<std::vector<double> > &, 
                         RanMars *, double, int,RDF *, bool);

// RDF
// take our RDF data, our value of nu, and the current step, and simulation #
void PrintRDF(RDF *, double, int, unsigned int);

// calculation of the FFT to get the power spectrum
void calculateFFT(RDF *);

int main () {
    
    // number of atoms in our system
    int nAtoms = 224;
    
    // equilibration steps to be run
    int equilibrationSteps = 1000000;

    // production steps to be run
    int productionSteps = 4000000;

    // our unit square simulation box
    Box box(1.0, 1.0);

    // double holding our current value of nu
    double nu;

    // double holding our current value of K
    double K;

    // seed for Marsaglia PRNG
    int seed = 28622498;
    
    // create the random number generator
    RanMars *thisRNG;
    thisRNG = new RanMars(seed);
   
    // initialize a pointer to our RDF class
    RDF *rdfData; 

    // initialize our atoms as a vector
    std::vector<Atom> atoms(nAtoms) ;
    
    // we make an NxN distance matrix
    // -- since we modify only one atom's position per step, it is inefficient to 
    //    calculate 0.5 * N^2 distances;
    //    We instead calculate ~N
    //    -- this will also rapidly speed up the RDF calculation and eliminate 
    //    redundancies, /especially/ when a move is rejected and we must compute the 
    //    /exact/ same RDF
    
    std::vector<std::vector<double> > distanceMatrix(nAtoms, std::vector<double>(nAtoms)) ;

    // our loop over all nu values, to do all the work iteratively
    for (unsigned int i = 0; i < (sizeof(nuValues) / sizeof(*nuValues)); i++ ) {
        // set our steps counter equal to zero
        steps = 0;
        // reset our 'accepted', 'total', and 'successRatio' variables for this simulation
        accepted = 0;
        total = 0;
        successRatio = 0;
        
        // initialize our atoms on the lattice for this iteration of the nu values
        // here, we also initialize the values stored in our distance matrix
        // so, to do so, we will need our box, which knows how to apply its PBCs
        InitializeCoordinates(atoms, box, distanceMatrix);

        /*
         * Look here for correct syntax of accessing the distance matrix;
         * note that the [i][j] syntax, i must be greater than j to get a meaningful distance
        // print the distance from atom 1 to atom 2;
        // this should be identically 1.0 / 14.0
        std::cout << "calculated distance atom 0 to atom 1: " << distanceMatrix[0][1] << std::endl;
        std::cout << "expected value:                       " << 1.0/14.0 << std::endl;

        std::cout << "calculated distance atom 0 to atom 14 " << distanceMatrix[0][14] << std::endl;
        double expectedDistance0_14 = sqrt( pow( (0.5 / 14), 2.0) + pow( (1.0/16.0), 2.0));
        std::cout << "expected value:                       " << expectedDistance0_14 << std::endl;
        */

        // set nu equal to nuValues[i]
        nu = nuValues[i];
        // set K equal to KValues[i]
        K = KValues[i];
        // compute the forbidden zone 'd0' based on the current value of nu in nuValues
        d0 = (1.0 / 14.0) * (1.0 - pow(2.0, nu - 8.0));
        // compute alpha for our given value of d0
        alpha = (1.0 / 14.0) - d0;
        
        // set the diameter of the atoms
        for (int j = 0; j < nAtoms; j++) {
            atoms[j].setDiameter(d0);    
        };
        // print the configuration; the pictures should all be atoms in the same position with varying size
        // this is the initial configuration set on the lattice
        PrintConfigurationData(atoms,nuValues[i],steps);
        
        
        // initialize our radial distribution function data
        // we will use 500 bins, with lowerBound 0.0 and upperBound 0.45 * boxLength
        rdfData = new RDF(500,0.0,K*d0, K);
        
        // do a number of Monte Carlo translations for 'equilibrationSteps' number of steps;
        // bool 'false' denoting production or not
        MonteCarloTranslate(atoms,box,distanceMatrix,thisRNG,nu,equilibrationSteps,rdfData,false);
       
        // do a number of Monte Carlo translations for 'productionSteps' number of steps;
        // bool 'true' denoting production run
        MonteCarloTranslate(atoms,box,distanceMatrix,thisRNG,nu,productionSteps,rdfData,true);

        // make a name for this data file unique to the simulation, specific to the value of nu
        PrintConfigurationData(atoms,nuValues[i],steps);
        
        // do the same with our rdfData
        // normalize the data first
        rdfData->normalizeRDF(productionSteps);
        PrintRDF(rdfData, nuValues[i], productionSteps, i);
        
        // calculate the free energy surface; F(r) = -kT ln (g(r));
        // we will then have F/kT = - ln (g(r))
        


        // make null, delete, and reseed the RNG for the next value of nu
        thisRNG = 0;
        delete thisRNG;
        // the numbers below were chosen completely arbitrarily; we just want different seeds
        thisRNG = new RanMars(seed+4 + 20*i);

        // and we do the same thing withour RDF data, after printing it
        rdfData = 0;
        delete rdfData;
        // we have a 'new' declaration in the next loop for rdfData; so we don't need it here
        
    };
    delete thisRNG;

    // exit the program
    return 0;


};

void InitializeCoordinates(std::vector<Atom> &atoms, Box &box, std::vector<std::vector<double> > &distances) {
    //int numAtoms = atoms.size();
    
    // the paper specifies we have a lattice 14x16; with spacing in the 
    // x-direction of 1/14, and spacing in the y direction of 1/16;
    // Every other row is offset by 0.5*(1/14) to create the offset lattice
    double initial_x_pos = 0.0 ; // we start the lattice at x = 0.0
    double initial_y_pos = 1.0 ; // we start the lattice at y = 1.0
    
    // we make the lattice atom by atom in the x-direction
    // after placing 14 particles in the row, we return to the initial x position
    // and place 14 more particles; we repeat this procedure 16 times.

    double x_increment = 1.0/14.0;
    double y_increment = 1.0/16.0;

    double x_offset = 0.5 * 1.0 / 14.0;

    double xpos = initial_x_pos;
    double ypos = initial_y_pos;

    int atomIndex = 0;

    // we specified that there are 16 rows and 14 columns..
    for (int i = 0; i < 16; i++) {
        // alternating rows are offset by x_offset; the first row is not
        if (i%2 == 1) {
            xpos = initial_x_pos + x_offset;
        } else {
            xpos = initial_x_pos;
        };
        // every row is y_increment distance below the other
        ypos = initial_y_pos - (i * y_increment);
        for (int j = 0; j < 14; j++) {
            atoms[atomIndex].setCoordinates(xpos,ypos);
            
            // increment the positions for the next atom to place in the lattice
            xpos += x_increment;
            
            // we don't increment y here; we are making the lattice row-by-row

            // increment the atomIndex
            atomIndex += 1;
        };
        
    }

    // populate distances for all m = [0,N-1], n = [i+1,N];
    // deliberately set all other indices to 0.0
    for (unsigned int m = 0; m < distances[0].size(); m++) {
        for (unsigned int n = 0; n < distances[0].size(); n++) {
            // if j <= i, set to zero; else compute the distance
            // -- we don't care about copy data that is redundant;
            // additionally, a particle is always distance 0.0 from itself
            if (n <= m) {
                distances[m][n] = 0.0;
            } else {
                distances[m][n] = box.computeDistance(atoms[m],atoms[n]);
            };
        };
    };

};

void MonteCarloTranslate(std::vector<Atom> &atoms, Box &thisBox, std::vector<std::vector<double> > &distanceMatrix,
                         RanMars *prng, double nu_, int steps_,RDF *rdf, bool production) {
    // create a vector with which to hold a list of coordinates, in the event that a move is rejected;
    // initialize it to zero
    std::vector<double> tempCoords = std::vector<double> (2, 0.0);

    // get the size of the atom list
    int nAtoms = atoms.size();

    // create a vector with which to hold a /vector/ of length N of old coordinates for some atom i between all other atoms j
    std::vector<double> pre_attempt_dist_ij = std::vector<double> (nAtoms, 0.0);
    
    std::vector<double> post_attempt_dist_ij = std::vector<double> (nAtoms, 0.0);
    // our number of atoms for indexing, since c++ is base 0
    int nAtomsIdx = nAtoms - 1;
    int idx;
    double dist;

    // flag denoting whether or not the move we are attempting is valid
    bool validMove;
    
    // if we are in production mode, initialize the RDF temporary data
    // (i.e., get the snapshot of the current configuration's rdf at this step
    // we will add the proper snapshot to gr at the end of the step (with modification if the move is accepted)
    if (production) {
        rdf->initialStep(distanceMatrix);
    };

    // for number of steps 'steps_', attempt a Monte Carlo translation move
    for (int i = 0; i < steps_; i++) {
        validMove = true;
        // randomly select an atom to translate; lrint rounds and casts as integer (from math.h)
        // we need nAtoms - 1, since prng->uniform goes [0,1) and max index of atoms is nAtoms-1
        idx = lrint (prng->uniform() * (nAtomsIdx));
        // store the atom's coordinates in the tempCoords vector, in case the move is rejected
        tempCoords = atoms[idx].coords;
        
        // ensure that pre_attempt_dist_ij and post_attempt_dist_ij are initialized to 0.0
        pre_attempt_dist_ij = std::vector<double> (nAtoms, 0.0);
        post_attempt_dist_ij = std::vector<double> (nAtoms, 0.0);

        // get the previous distances b/w atom idx and all other atoms j;
        for (int m = 0; m < nAtoms; m++) {
            if (m < idx) {
                pre_attempt_dist_ij[m] = distanceMatrix[m][idx];
            } else if (m > idx) {
                pre_attempt_dist_ij[m] = distanceMatrix[idx][m];
            } else {
                pre_attempt_dist_ij[m] = 0.0;
            };
        };

        
        // randomly displace the atom
        thisBox.Translate(atoms[idx], alpha, prng);
       
        // increment the total number of steps taken in simulation
        total += 1.0;

        //thisBox.Translate(atoms[idx],alpha);
        
        // accept or reject the move
        // essentially, reject if the distance is less than one diameter
        for (int j = 0; j < nAtoms; j++) {
            if (j != idx) {
                dist = thisBox.computeDistance(atoms[idx],atoms[j]);
                post_attempt_dist_ij[j] = dist;
                if (dist < atoms[idx].diameter) {
                    validMove = false;
                    // also, we can exit the for loop
                    break;
                };
            } else {
                post_attempt_dist_ij[j] = 0.0;
            };
        };

        // if we accept the move, do stuff;
        // otherwise, reset the coordinates of the attempted move
        if (validMove) {
            
            //update distance matrix; the only values perturbed are those concerning atom idx
            for (int ii = 0; ii < nAtoms; ii++) {
                if (ii < idx) {
                    distanceMatrix[ii][idx] = post_attempt_dist_ij[ii];
                } else if (ii > idx) {
                    distanceMatrix[idx][ii] = post_attempt_dist_ij[ii];
                } else {
                    distanceMatrix[idx][ii] = 0.0;
                };

            };
            accepted += 1.0;

        } else {
            // return to its original position
            atoms[idx].coords = tempCoords;
        };
        successRatio = accepted / total;

        if (i%100000 == 0) {
            printf("At step %d; production mode %d\n",i,production);
            printf("Success ratio currently %f\n", successRatio);
        };

        if (production) {
            if (validMove) {
                // we need to update our snapshot of the radial distribution function
                // pass in the pre-attempt vector of distances, post-attempt, and idx of perturbed atom
                rdf->updateSnapshot(pre_attempt_dist_ij, post_attempt_dist_ij, idx);
            };
            //rdf->calculateRDF(atoms,thisBox);
            rdf->calculateRDF();
            steps += 1;
            //if (i%10000 == 0 && i!=0) {
                //PrintConfigurationData(atoms,nu_,i);
            //};

        };
    };
};

void PrintConfigurationData(std::vector<Atom> &atoms, double nu_, int numSteps) {

    std::string xyzFileName = "";
    // we need to make a file name unique to this nu value and step number;
    // flush the stream, in case other stuff was in there
    std::ostringstream stringStream;
    stringStream.flush();
    stringStream.str("");
    stringStream << "Nu_" << nu_ << "_step_" << numSteps << ".xyz";
    xyzFileName = stringStream.str();
    // convert xyzFileName from std::string to const char
    std::ofstream xyzFile(xyzFileName.c_str(), std::ios::out);

    double hardSphereRadius = 0.5 * atoms[0].diameter;
    // the xyz file will contain the x, y, and radius, in that order
    int numAtoms = atoms.size();
    xyzFile << "# Hard sphere simulation of " << numAtoms << " atoms with radius " << hardSphereRadius << std::endl;
    xyzFile << "#\n#\n" << "# format: xpos ypos radius" << std::endl;

    for (unsigned int i = 0; i < atoms.size(); i++) {
        xyzFile << atoms[i].coords[0] << " " << atoms[i].coords[1] << " " << hardSphereRadius << std::endl;
    };

};

// takes the rdf, the value of nu, and the number of steps that we have taken;
// and also, the simulation iteration (for when we do duplicate nu/K values
// to reproduce the original paper
void PrintRDF(RDF *rdfData,  double nu_, int numSteps, unsigned int simIteration) {

    std::string rdfFileName = "";
    std::ostringstream stringStream;
    stringStream.flush();
    stringStream.str("");
    stringStream << "Nu_" << nu_ << "_step_" << numSteps << "sim" << simIteration << ".rdf";

    rdfFileName = stringStream.str();
    std::ofstream rdfFile(rdfFileName.c_str(), std::ios::out);

    // the number of bins is held in rdfData's member 'M'
    // alternatively, we might get the .size() of rdfData->gr
    int nBins = rdfData->M;

    // initialize the bin center to be lowerBound + 0.5 * dr;
    double binCenter = (rdfData->lowerBound) + (0.5 * rdfData->dr);
    
    // we must normalize this by the forbidden radius specified by this nu value;
    binCenter /= rdfData->d0;
    // nondimensionalize the binCenter - i.e., divide it by lowerBound (d0) and increment it 
    // by dr/d0
    for (int i = 0; i<nBins; i++) {
        rdfFile << binCenter << " " << rdfData->gr[i] << std::endl;

        // increment the bin by dr
        binCenter += (rdfData->dr / rdfData->d0);
    };
};




