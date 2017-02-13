#include "simulation.h"
#include <iostream>
#include <math.h>
#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>

// declare our constructor
Simulation::Simulation(int _numberOfAtoms, double _density, double _dt, double _mass,
                   double _sigma, double _epsilon, double _Tstar, int _seed, 
                   std::string _name) :
            numberOfAtoms(_numberOfAtoms), density(_density), dt(_dt), mass(_mass),
            sigma(_sigma), epsilon(_epsilon), Tstar(_Tstar), name(_name)
{
// and stuff to do here
    prng = new RanMars(_seed);
    box = new Box(numberOfAtoms, density, sigma);
    step = 0;
    
    rcut = 3.0 * sigma;
    sig3 = sigma * sigma * sigma;
    sig6 = sig3 * sig3;
    sig12 = sig6 * sig6;

    rcut3 = rcut * rcut * rcut;
    rcut9 = rcut3 * rcut3 * rcut3;
    tail = 0.0;

    dt2m = dt / (2.0 * mass);
    dt2 = dt / 2.0;
    
    deltaR2 = std::vector<double> (numberOfAtoms, 0.0);
};


void Simulation::initializeAtomsPositions() {
    // initialize the atoms on the lattice in the box;

    // initialize the vector of atoms
    atoms = std::vector<Atom> (numberOfAtoms);

    // and tell box to place the atoms in the box & compute distances subject to the geometry
    box->initializeAtoms(atoms);


}; 

void Simulation::initializeAtomsVelocities() {

    for (int i = 0; i < (int) atoms.size(); i++) {
        // we need to sample from the proper distribution here corresponding to our T*...


    };
    
    
    
    // calculate the total momentum
    //calculateTotalMomentum();
    //
    //scale the momentum s.t. that net momentum within the box is zero
    //scaleMomentum();

};
// our run() method, for nsteps;
// print xyz every printXYZ steps,
// and a bool parameter denoting whether or not we are running production
void Simulation::run(int nsteps, int printXYZ, bool production) {
    
    std::ostringstream stringStream;
    stringStream.flush();
    stringStream.str("");
    stringStream << name << "_prop" << ".dat";
    std::string simData = stringStream.str();
    std::ofstream simDataFile(simData.c_str(), std::ios::out);
    
    simDataFile << "# format: step  potentialEnergy kineticEnergy totalEnergy" << std::endl;
    total = 0.0;

    int nAtomsIdx = atoms.size() - 1;
    int nAtoms = atoms.size();

    //std::cout << "In run() method; about to call ComputeTotalEnergy() " << std::endl;
    // get the current potential energy of the system; also initializes the tail correction
    ComputeTotalPotentialEnergy();


    // we are initializing the simulation; get the current force on each atom
    ComputeTotalForce();

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
        for (int ijk = 0; ijk < (int) atoms.size(); ijk++) {
            forcesToPrint = atoms[ijk].getTotalForce();
            initForcesFile << ijk << " " << forcesToPrint[0] << " " << forcesToPrint[1] << " " 
                << forcesToPrint[2] << std::endl;
        };
    };

    for (int i = 0; i < nsteps; i++) {
        
        /* Velocity-Verlet Integrator
         *
         * --update velocities by half step, according to current forces (accelerations)
         * --update positions, according to current velocities & accelerations
         * --zero the forces
         * --compute the forces present in the system
         * --update the velocities by a half step
         * ----and do this 'nsteps' number of times
         *
         */

        // half-step velocities update
        velocityHalfstep();

        // full step positions update
        // --note: apply PBC's within this for each atom as it is translated!
        positionStep();

        // zero the fs vector for each atom
        zeroForces();

        // re-calculate the forces
        calculateForces();

        // update the velocities by a halfstep
        velocityHalfstep();

        // take the average of the deltaR2 vector right now
        calculateAvgDeltaR2();

        // if we are at the end of printEvery during a production cycle, do stuff:
        if ( ( (i % printXYZ ) == 0 ) and (production) ) {
            std::cout << "At step " << i << " of production." << std::endl;
            printConfig(name,i);           
        };

        // output to simData the following quantities each step:
        // step   PE    KE   (totalE)     
        double totalPE_n, totalE_n;
        // we print our properties out normalized s.t. they are per-atom
        totalPE_n = totalPE / ((double) atoms.size());
        totalE_n = totalPE_n + totalKE_n;

        simDataFile << i  << " " << totalPE_n << " " << totalKE_n << " " << (totalE_n) << std::endl;
    };
};

void Simulation::velocityHalfStep() {
    // here, we modify the velocities by a halfstep according to the velocity-verlet integration scheme
    std::vector<double> force, dv;

    // NOTE: at this point, each atom should know all pairwise forces acting on it;
    // they just need to be integrated s.t. the velocity is modified accordingly
    for (int i = 0; i < numberOfAtoms; i++) {
        force = atoms[i].getForces();
        dv[0] = dt2m * force[0];
        dv[1] = dt2m * force[1];
        dv[2] = dt2m * force[2];

        atoms[i].addVelocities(dv);
    };
};

void Simulation::positionStep() {
    std::vector<double> vel, dx;

    double dr2; // deltaR^2
    // NOTE: at this point, each atom should know its velocity;
    // we simply need to integrate it s.t. the position is modified accordingly
    // ALSO, here, we might keep track of the < ( delta r)^2 > quantity;
    // additionally, we would enforce PBC in this function /after/ modifying the position

    for (int i = 0; i < numberOfAtoms; i++) {
        vel = atoms[i].getVelocities();
        dx[0] = dt * vel[0];
        dx[1] = dt * vel[1];
        dx[2] = dt * vel[2];

        dr2 = (dx[0] * dx[0]) + (dx[1] * dx[1]) + (dx[2] * dx[2]);

        deltaR2[i] += dr2;

        atoms[i].addPositions(dx);

        box->enforcePBC(atoms[i]);
    };
};

void Simulation::zeroForces() {
    std::vector<double> zeros = std::vector<double> (3, 0.0);

    for (int i = 0; i < numberOfAtoms; i++) {
        atoms[i].setForces(zeros);
    };
};

// update fs vector for each atom
void Simulation::calculateForces() {
    for (int i = 0; i < (numberOfAtoms - 1); i++) {
        for (int j = i+1; j < numberOfAtoms; j++) {
            LJForce(i,j);
        };
    };
};

void Simulation::ComputeTotalPotentialEnergy() {
    // zero our totalPE variable
    totalPE = 0.0;
    for (int i = 0; i < (numberOfAtoms - 1); i++) {
        for (int j = i+1; j < (int) numberOfAtoms; j++) {
            totalPE += LJPotential(i,j);
        };
    };

    // and our tail correction, since we truncate the potential abruptly
    // note that 'density' here is rhoStar - so, divide by sig3
    tail = ( (8.0 / 3.0) * M_PI * density / (sig3));
    tail *= ( ( (1.0 / 3.0) * (sig9 / rcut9)) - (sig3 / rcut3));
    totalPE += tail;
};

void Simulation::ComputeTotalForce() {
    // populate forces for all atoms in the simulation 
    // (note: we print this immediately after initializing the lattice..)

    //std::cout << "in ComputeTotalForce!" << std::endl;
    for (int i = 0; i < (int) (numberOfAtoms - 1); i++) {
        for (int j = i + 1; j < (int) numberOfAtoms ; j++) {
            LJForce(i,j);
        };
    };
};

double Simulation::LJPotential(int ii,int jj) {

    double rij = box->computeDistance(atoms[ii],atoms[jj])
    
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
    
    double rij = box->computeDistance(atoms[ii],atoms[jj]);

    std::vector<double> force_vec = std::vector<double> (3, 0.0);
    std::vector<double> force_vec_neg = force_vec;
    // we will need the unit vector
    
    if (rij < rcut) {
        std::vector<double> coords1 = atoms[i].getCoordinates();
        std::vector<double> coords2 = atoms[j].getCoordinates();
        std::vector<double> rij_unit_vec = std::vector<double> (3, 0.0);

        double rij3 = rij * rij * rij;
        double rij6 = rij3 * rij3;
        double rij12 = rij6 * rij6;
        
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

        atoms[i].addForce(force_vec);
        atoms[j].addForce(force_vec_neg);

    } else {
        // do nothing; simply return
        return;
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


