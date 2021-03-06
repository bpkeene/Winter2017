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
                   double _sigma, double _epsilon, double _Tstar, int _seed, 
                   std::string _name) :
            numberOfAtoms(_numberOfAtoms), density(_density), 
            sigma(_sigma), epsilon(_epsilon), Tstar(_Tstar), name(_name)
{
// and stuff to do here
    prng = new RanMars(_seed);
    box = new Box(numberOfAtoms, density, sigma);
    step = 0;
    distances = std::vector< std::vector<double> > (numberOfAtoms, std::vector<double> (numberOfAtoms));
    ijDistance = std::vector<double> (numberOfAtoms, 0.0);
    ijDistance_old = std::vector<double> (numberOfAtoms, 0.0);
    alpha = 1.0;
    rcut = 3.0 * sigma;
    sig3 = sigma * sigma * sigma;
    sig6 = sig3 * sig3;
    sig12 = sig6 * sig6;

    rcut3 = rcut * rcut * rcut;
    rcut9 = rcut3 * rcut3 * rcut3;
    tail = 0.0;

};


void Simulation::initializeAtoms() {
    // initialize the atoms on the lattice in the box;
    // what should be the inter-atomic distance?

    // initialize the vector of atoms
    atoms = std::vector<Atom> (numberOfAtoms);

    // and tell box to place the atoms in the box & compute distances subject to the geometry
    box->initializeAtoms(atoms,distances);
    //std::cout << "returned from box->initializeAtoms" << std::endl;
}; 

// our run() method, for nsteps;
// print xyz every printXYZ steps,
// and a bool parameter denoting whether or not we are running production
void Simulation::run(int nsteps, int printXYZ, bool production) {
    double rho = (double) numberOfAtoms * sig3 / (box->getVolume());
    double rho2 = rho * rho;
    // calculate the pressure tail; this will be a static quantity throughout the simulation
    ptail = (16.0 / 3.0) * M_PI * rho2 * ( (2.0/3.0) * (sig3*sig6 /rcut9) - (sig3/rcut3));
    ptail_n = ptail / ( (double) numberOfAtoms);

    // calculate the kinetic contribution to the pressure
    pkinetic = rho * epsilon * Tstar;
    pkinetic_n = pkinetic / ( (double) numberOfAtoms);
    
    totalKE = 3.0 / 2.0 *  ( (double) numberOfAtoms) * epsilon * Tstar;
    totalKE_n = 3.0 / 2.0 * epsilon * Tstar;
    
    // if we are doing a production run, open a file in which we can write data throughout the simulation
    // we'll also write a configuration file every 'printXYZ' number of steps
    
    std::ostringstream stringStream;
    stringStream.flush();
    stringStream.str("");
    stringStream << name << "_prop" << ".dat";
    std::string simData = stringStream.str();
    std::ofstream simDataFile(simData.c_str(), std::ios::out);
    if (production) {
        simDataFile << "# format: step    potentialEnergy   kineticEnergy    totalEnergy   pKinetic   pVirial   pTot" << std::endl;
    };
    // for a given number of steps..
    accepted = 0.0;
    total = 0.0;

    // and we'll do block averages as well for sampling for good alpha
    double acceptedInterval = 0.0;
    double totalInterval = 0.0;

    int idx; // our randomly selected index - the atom which we will move
    
    double pre_move_pot = 0.0;
    double post_move_pot = 0.0;
  
    // every 5% of moves, change alpha by some factor if local acceptance is not 40 < x < 60 %
    int changeAlphaEvery = (int) floor(0.05 * nsteps);

    int nAtomsIdx = atoms.size() - 1;
    int nAtoms = atoms.size();

    //std::cout << "In run() method; about to call ComputeTotalEnergy() " << std::endl;
    // get the current potential energy of the system; also initializes the tail correction
    ComputeTotalEnergy();

    if (! (production) ) {
        for (int qq = 0; qq < atoms.size(); qq++ ) {
            atoms[qq].initializeForces(nAtoms);
        };
    };

    //std::cout << "Called ComputeTotalEnergy(); got a value of " << totalPE << std::endl;

    ComputeTotalForce();

    //std::cout << "Called ComputeTotalForce() successfully" << std::endl;
    // and populate the totalForce vector within each atom
    // -- this amounts to adding up N vectors of size Nx3.. but we do it once, 
    // and then every operation is O(N) afterwards
    for (int q = 0; q < atoms.size(); q++) {
        atoms[q].computeTotalForce();
    };

    // if we are doing our equilibration, we can now write the forces data to file
    // -- note that we expect the net force on each atom to be zero
    if (!(production)) {
        std::ostringstream stringStream2;
        stringStream2.flush();
        stringStream2.str("");
        stringStream2 << name << "_init_forces.dat";
        std::string initForces = stringStream2.str();
        std::ofstream initForcesFile(initForces.c_str(), std::ios::out);
        initForcesFile << "# format: atomIndex Fx Fy Fz" << std::endl;
        initForcesFile << atoms.size() << "\n" << std::endl;
        std::vector<double> forcesToPrint;
        for (int ijk = 0; ijk < atoms.size(); ijk++) {
            forcesToPrint = atoms[ijk].getTotalForce();
            initForcesFile << ijk << " " << forcesToPrint[0] << " " << forcesToPrint[1] << " " 
                << forcesToPrint[2] << std::endl;
        };
    };

    // just call this method anyways..
    ComputePressure();

    //std::cout << "About to begin a run!" << std::endl;
    for (int i = 0; i < nsteps; i++) {
        //std::cout << "on step " << i << std::endl;
        // randomly select an atom to move
        idx = lrint(prng->uniform() * nAtomsIdx);
        
        // save its coordinates and forces vector
        atoms[idx].saveCoordinates();
        atoms[idx].saveForces();
        
        // all other atoms should likewise save their forces w.r.t. atom idx
        for (int jj = 0; jj < atoms.size(); jj++) {
            if (jj != idx) {
                atoms[jj].saveForce(idx);
            };
        };

        // extract the temporary distance vector;
        extractTempDistanceVector(idx);

        // compute the energy of its current coordinates
        pre_move_pot = ComputeAtomPotential(idx);
        
        // save the temporary distance vector
        saveTempDistanceVector();
        
        // translate the atom, applying PBC
        box->Translate(atoms[idx], alpha, prng); 

        // ComputeAtomPotential assumes our distances matrix is up-to-date...
        // alter the values in our distances matrix 
        // and put these values in to ijDistances
        
        // update the distances matrix 
        updateDistancesMatrix(idx);

        // update the temporary distance vector
        extractTempDistanceVector(idx);

        // compute delta E
        post_move_pot = ComputeAtomPotential(idx);

        // accept or reject the move
        acceptMove = acceptOrReject(pre_move_pot, post_move_pot);

        if (acceptMove) {
            acceptedInterval += 1.0;
            accepted += 1.0;

            totalPE += (post_move_pot - pre_move_pot);
            for (int kk = 0; kk < atoms.size(); kk++) {
                if (kk != idx) {
                    // update the forces on each atom to be consistent with the new position
                    LJForce(kk,idx);
                    atoms[kk].updateTotalForce(idx);
                };
            };

            ComputeDeltaPressure(idx);

        } else {
            
            // set the coordinates of this atom back to its old coordinates
            atoms[idx].resetCoordinates();
            // reset the distances in our distance matrix
            // --need to know what the reference atom is!
            resetDistancesMatrix(idx);
            atoms[idx].resetForces();
            
            for (int jjj = 0; jjj < atoms.size(); jjj++) {
                if (jjj != idx) {
                    atoms[jjj].resetForce();
                };
            };
        };

        // update the total move counters, for this interval and for the simulation
        totalInterval += 1.0;
        total += 1.0;


        // if we are at the end of a block, i!=0, and doing equilibration, check alpha
        if ( ( (i % changeAlphaEvery) == 0) and (i != 0) and (!(production))) {
            double successRatioInterval = acceptedInterval / totalInterval;
            if (successRatioInterval >= 0.55) {
                // we can increase alpha a bit to take larger steps
                alpha += (prng->uniform() * prng->uniform()) * sigma;
            }; 
            if (successRatioInterval <= 0.45) {
                // we should take smaller steps.. note that uniform is [0,1.) so we don't need fabs()
                alpha -= (prng->uniform() * prng->uniform())*alpha;
            };
            // sanity check: print the equilibration step # (absolute) and the interval acceptance ratio
            std::cout << "At step " << i << " of equilibration " << "with acceptance ratio " << successRatioInterval << std::endl;
            // reset the counters for the next interval
            acceptedInterval = 0.0;
            totalInterval = 0.0;
        };

        // if we are at the end of printEvery during a production cycle, do stuff:
        if ( ( (i % printXYZ ) == 0 ) and (production) ) {
            std::cout << "At step " << i << " of production." << std::endl;
            printConfig(name,i);           
        };

        if (production) {
        // output to simData the following quantities each step:
        // step   PE    KE   (totalE)   pKinetic     PVirial   PTot   
            double totalPE_n, totalE_n;
        // we print our properties out normalized s.t. they are per-atom
            totalPE_n = totalPE / ((double) atoms.size());
            totalE_n = totalPE_n + totalKE_n;



            simDataFile << i  << " " << totalPE_n << " " << totalKE_n << " " << (totalE_n) << 
            " " << pkinetic << " " << pvirial << " " << (pkinetic + pvirial + ptail) << std::endl;
        };
    };
};


bool Simulation::acceptOrReject(double pre, double post) {
    double deltaE = post - pre;
    if (deltaE < 0.0) {
        return true;
    } else {
        double probability = exp(-deltaE/(epsilon*Tstar));
        double rn = prng->uniform();
        if (probability > rn) {
            return true;
        } else {
            return false;
        };
    };
};

void Simulation::saveTempDistanceVector() {
    ijDistance_old = ijDistance;
};


void Simulation::ComputePressure() {
    
    double volume = box->getVolume();

    double virialSum = 0.0;

    // we get the actual vector r_ij;
    std::vector<double> coords1;
    std::vector<double> coords2;

    // the magnitude of r_ij, as stored in our distances matrix
    double rij;
    double rijx, rijy, rijz;
    std::vector<double> forces_ij;

    double fs;

    // do the sum
    for (int i = 0; i < (atoms.size() - 1); i++) {
        coords1 = atoms[i].getCoordinates();
        for (int j = i+1; j < atoms.size(); j++) {
            fs = 0.0;
            coords2 = atoms[j].getCoordinates();
            rij = distances[i][j];
            rijx = (coords2[0] - coords1[0]) / rij;
            rijy = (coords2[1] - coords1[1]) / rij;
            rijz = (coords2[2] - coords1[2]) / rij;
            forces_ij = atoms[i].getForces(j);

            fs = (forces_ij[0] * rijx) + (forces_ij[1] * rijy) +
                 (forces_ij[2] * rijz);

            virialSum += fs;
        };
    };
    double N = (double) atoms.size();
    pvirial = virialSum * 2.0 / (3.0 * volume );
};

// after initially calculating the pressure in the simulation, we just need delta P's
// from move to move
void Simulation::ComputeDeltaPressure(int idx) {
    double volume = box->getVolume();

    double deltaVirial = 0.0;
    double oldPressure = 0.0; 
    double newPressure = 0.0;

    double sumDeltaPs = 0.0;
    // we get the actual vetor_r_ij;
    std::vector<double> coords1;
    std::vector<double> coords2;

    std::vector<double> coords1_o;

    // the magnitude of r_ij, as stored in our distances matrix
    double rij, rij2;
    double rijx, rijy, rijz;
    double rijx2, rijy2, rijz2;

    // get this atom's forces vector from before it was perturbed
    std::vector< std::vector<double> > oldForces = atoms[idx].getOldForces();
   
    // we can assume the following:
    std::vector< std::vector<double> > newForces = atoms[idx].getForcesMatrix();
    // all forces on every atom are up to date, as are the distances
    // also, the distances are up to date
    coords1 = atoms[idx].getCoordinates();
    coords1_o = atoms[idx].getOldCoordinates();  

    int ii, jj;
    for (int i = 0; i < (atoms.size()); i++) {
        coords2 = atoms[i].getCoordinates();
        if (i < idx) {
            ii = i;
            jj = idx;
        } else if (i > idx) {
            ii = idx;
            jj = i;
        } else {
            continue;
        };
        rij = distances[ii][jj];
        rij2 = ijDistance_old[i];
        
        rijx = (coords2[0] - coords1[0]) / rij;
        rijy = (coords2[1] - coords1[1]) / rij;
        rijz = (coords2[2] - coords1[2]) / rij;
    
        rijx2 = (coords2[0] - coords1_o[0]) / rij2;
        rijy2 = (coords2[1] - coords1_o[1]) / rij2;
        rijz2 = (coords2[2] - coords1_o[2]) / rij2;

        oldPressure = (oldForces[i][0]*rijx2) + (oldForces[i][1]*rijy2) + (oldForces[i][2]*rijz2);
        newPressure = (newForces[i][0]*rijx ) + (newForces[i][1]*rijy ) + (newForces[i][2]*rijz );
        sumDeltaPs += (newPressure - oldPressure);
    };

    double N = (double) atoms.size();
    deltaVirial = sumDeltaPs * 2.0 / ( 3.0 * volume );

    pvirial += deltaVirial;

};
void Simulation::extractTempDistanceVector(int idx) {
    int ii, jj;
    for (unsigned int i = 0; i < atoms.size(); i++) {
        if (i < idx) {
            ii = i;
            jj = idx;
        } else if (i > idx) {
            ii = idx;
            jj = i;
        } else {
            ijDistance[i] = 0.0; // we are getting the distance of an atom from itself
            continue;
        };
        
        ijDistance[i] = distances[ii][jj];
    };
};

void Simulation::updateDistancesMatrix(int idx) {
    // atom at idx was perturbed; for all atoms j, calculated distance idx-j
    int ii, jj;
    for (int i = 0; i < distances[0].size(); i++) {
        if (i < idx) {
            ii = i;
            jj = idx;
        } else if (i > idx) {
            ii = idx;
            jj = i;
        } else {
            distances[i][idx] = 0.0;
            continue;
        };
        distances[ii][jj] = box->computeDistance(atoms[ii],atoms[jj]);
    };
};

void Simulation::resetDistancesMatrix(int idx) {
    // atom at idx was perturbed, but move was rejected;
    int ii, jj;
    for (int i = 0; i < distances[0].size(); i++) {
        if (i < idx) {
            ii = i;
            jj = idx;
        } else if (i > idx) {
            ii = idx;
            jj = i;
        } else { 
            distances[i][idx] = 0.0;
            continue;
        };
        distances[ii][jj] = ijDistance_old[i];
    };
};

void Simulation::ComputeTotalEnergy() {
    // zero our totalPE variable
    totalPE = 0.0;
    for (int i = 0; i < (distances[0].size() - 1); i++) {
        for (int j = i+1; j < distances[0].size(); j++) {
            totalPE += LJPotential(i,j);
        };
    };

    // and our tail correction, since we truncate the potential abruptly
    // note that 'density' here is rhoStar - so, divide by sig3
    tail = ( (8.0 / 3.0) * M_PI * density / (sig3));
    tail *= ( ( (1.0 / 3.0) * (1.0 / rcut9)) - (1.0 / rcut3));
    totalPE += tail;
};

void Simulation::ComputeTotalForce() {
    // populate forces for all atoms in the simulation 
    // (note: we print this immediately after initializing the lattice..)

    //std::cout << "in ComputeTotalForce!" << std::endl;
    for (int i = 0; i < (atoms.size() - 1); i++) {
        for (int j = i + 1; j < atoms.size() ; j++) {
            //std::cout << "Computing LJForce between " << i << " and " << j << std::endl;
            LJForce(i,j);
        };
    };
};

double Simulation::ComputeAtomPotential(int idx) {
    double pe = 0.0;
    for (int i = 0; i < atoms.size(); i++) {
        pe += LJPotential(i,idx);
    };
    return pe;
};


double Simulation::LJPotential(int ii,int jj) {
    int i, j;
    // set i to the lesser index, so that we access the proper 
    // data in distances matrix
    if (ii < jj) {
        i = ii;
        j = jj;
    } else if (ii > jj) {
        i = jj;
        j = ii;
    } else {
        return 0.0;
    };

    double rij = distances[i][j];
    
    // note that we mixed up ii, jj;
    // but energy is invariant to this (not a directional quantity);
    // really only because this quantity is a function of /distance/, and not directional
    if (rij < rcut) {
        // we have global values for sigma, eps, sig6, sig12...

        double ijEnergy = 0.0;

        double rijsq = rij * rij;
        double rij6 = rijsq * rijsq * rijsq;
        double rij12 = rij6 * rij6;

        ijEnergy = 4.0 * epsilon * ( (sig12/rij12) - (sig6/rij6) ); 

        return ijEnergy;
    } else {
        return 0.0;
    };

};

void Simulation::LJForce(int ii, int jj) {
    int i, j;
    // set i to the lesser index, so that we access the proper
    // data in distances matrix
    // initialize some vector to return to 0.0's for Fx, Fy, Fz
    std::vector<double> ijForce = std::vector<double> (3, 0.0);
    if (ii < jj) {
        i = ii;
        j = jj;
    } else if (ii > jj) {
        i = jj;
        j = ii;
    } else {
        return;    
    };
   
    //std::cout << "in LJForce, about to get distances[i][j]" << std::endl;

    // first, get the distance between the two atoms..
    double rij = distances[i][j];

    //std::cout << "value of rij: " << rij << std::endl;

    std::vector<double> force_vec = std::vector<double> (3, 0.0);
    std::vector<double> force_vec_neg = force_vec;
    // we will need the unit vector
    
    if (rij < rcut) {
        std::vector<double> coords1 = atoms[i].getCoordinates();
        std::vector<double> coords2 = atoms[j].getCoordinates();
        std::vector<double> rij_unit_vec = std::vector<double> (3, 0.0);

        //std::cout << "rij is less than rcut! we're going to compute this. " << std::endl;
        double rij3 = rij * rij * rij;
        double rij6 = rij3 * rij3;
        double rij12 = rij6 * rij6;
        
        //std::cout << 
        // we will calculate the magnitude of the force and dot it with rij_unit_vec
        // to get its directionality
        rij_unit_vec[0] = (coords2[0] - coords1[0]) / rij;
        rij_unit_vec[1] = (coords2[1] - coords1[1]) / rij;
        rij_unit_vec[2] = (coords2[2] - coords1[2]) / rij;

        // magnitude of our force from atom i to atom j
        double ijForce = (24.0 * epsilon / rij) * ( 2.0 * (sig12 / rij12) - (sig6 / rij6) ) ;

        force_vec[0] = ijForce * rij_unit_vec[0];
        force_vec[1] = ijForce * rij_unit_vec[1];
        force_vec[2] = ijForce * rij_unit_vec[2];

        force_vec_neg[0] = - force_vec[0];
        force_vec_neg[1] = - force_vec[1];
        force_vec_neg[2] = - force_vec[2];


        // we now have force;
        // we should have saved the force for each atom at some point, so just set it in each vector
        atoms[i].setForce(j,force_vec);
        atoms[j].setForce(i,force_vec_neg);

        // do we need to return this? idk..
    } else {
        // force_vec
        atoms[i].setForce(j,force_vec);
        atoms[j].setForce(i,force_vec);
    };
    return;
};

void Simulation::printConfig(std::string name, int step) {
    std::ostringstream stringStream;
    stringStream.flush();
    stringStream.str("");
    stringStream << name << "_step_" << step << ".xyz";
    std::string xyzFileName = stringStream.str();
    std::ofstream xyzFile(xyzFileName.c_str(), std::ios::out);
    
    // write the number of atoms, and then a comment line denoting the format
    xyzFile << atoms.size() << std::endl;
    xyzFile << "# format: xpos ypos zpos" << std::endl;
   
    // now, write atom name, xpos, ypos, zpos (newline) to file
    std::vector<double> theseCoords = std::vector<double> (3,0.0);
    
    for (unsigned int i = 0; i < atoms.size(); i++) {
        theseCoords = atoms[i].getCoordinates();
        xyzFile << "atom" << i << " " << theseCoords[0] << " " << theseCoords[1] << " " << theseCoords[2] << std::endl;
    };

};


