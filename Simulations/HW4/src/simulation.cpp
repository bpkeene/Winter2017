#include "simulation.h"
#include <iostream>
#include <math.h>
#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <stdlib.h>
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
    
    rcut = 3.0;
    sig2 = sigma * sigma;
    sig3 = sigma * sigma * sigma;
    sig6 = sig3 * sig3;
    sig9 = sig3 * sig3 * sig3;
    sig12 = sig6 * sig6;

    rcut3 = rcut * rcut * rcut;
    rcut9 = rcut3 * rcut3 * rcut3;
    tail = 0.0;
    double doubleCastMass = (double) mass;
    dt2m = dt / (2.0 * doubleCastMass);
    dt2 = dt / 2.0;
    
    deltaR2 = std::vector<double> (numberOfAtoms, 0.0);
}


void Simulation::initializeAtomsPositions() {
    // initialize the atoms on the lattice in the box;

    // initialize the vector of atoms
    atoms = std::vector<Atom> (numberOfAtoms);

    // and tell box to place the atoms in the box & compute distances subject to the geometry
    box->initializeAtoms(atoms,mass);
} 

void Simulation::initializeAtomsVelocities() {

    // these will be obtained from our prng
    double vx, vy, vz;

    // our cumulative momenta in the x, y, z directions
    double px, py, pz;

    // to normalize the momenta
    double nAtoms = (double) atoms.size();

    // cast the mass as a temporary variable, type double
    double massTemporary = (double) mass;

    // our BoxMullerConstant, to be used in generating the velocities corresponding to our T*
    //double BoxMullerConstant = sqrt(epsilon * Tstar / massTemporary);
    
    px = py = pz = 0.0;

    std::vector<double> theseVels;
    for (int i = 0; i < (int) atoms.size(); i++) {
        Atom &a = atoms[i];
        // we need to sample from the proper distribution here corresponding to our T*...
        vx = prng->gaussian() * Tstar;
        vy = prng->gaussian() * Tstar;
        vz = prng->gaussian() * Tstar;
        a.setVelocities(vx,vy,vz);
        theseVels = a.getVelocities();
        px += (theseVels[0] );
        py += (theseVels[1] );
        pz += (theseVels[2] );
    };
 
    px /= nAtoms;
    py /= nAtoms;
    pz /= nAtoms;


    double sumvSqr = 0.0;

    for (int k = 0; k < (int) atoms.size(); k++) {
        Atom &b = atoms[k];
        theseVels = b.getVelocities();
        theseVels[0] -= (px );
        theseVels[1] -= (py );
        theseVels[2] -= (pz );
        b.setVelocities(theseVels[0], theseVels[1], theseVels[2]);
        sumvSqr += 2.0 * (massTemporary * 0.5 * ( (theseVels[0] * theseVels[0]) + 
                                                  (theseVels[1] * theseVels[1]) + 
                                                  (theseVels[2] * theseVels[2])));
    };

    sumvSqr /= nAtoms;
    
    for (int m = 0; m < (int) atoms.size(); m++) {
        Atom &d = atoms[m];
        theseVels = d.getVelocities();
        theseVels[0] *= sqrt(3.0 * Tstar / sumvSqr);
        theseVels[1] *= sqrt(3.0 * Tstar / sumvSqr);
        theseVels[2] *= sqrt(3.0 * Tstar / sumvSqr);
        d.setVelocities(theseVels[0], theseVels[1], theseVels[2]);
    }


    // and finally, verify that our total momentum is now zero in the three directions
    px = py = pz = 0.0;
    for (int ijk = 0; ijk < (int) atoms.size(); ijk++) {
        Atom &c = atoms[ijk];
        theseVels = c.getVelocities();
        px += theseVels[0] * massTemporary;
        py += theseVels[1] * massTemporary;
        pz += theseVels[2] * massTemporary;
    };

    std::cout << "Momenta px py pz: " << px << " " << py << " " << pz << std::endl;

}


// our run() method, for nsteps;
// print xyz every printXYZ steps,
// and a bool parameter denoting whether or not we are running production
void Simulation::run(int equilSteps,int prodSteps, int printXYZ, bool production) {
    
    //std::cout << "made it here" << std::endl;
    std::ostringstream stringStream;
    stringStream.flush();
    stringStream.str("");
    stringStream << name << "_prop" << ".dat";
    std::string simData = stringStream.str();
    std::ofstream simDataFile(simData.c_str(), std::ios::out);
    
    simDataFile << "# format: step  potentialEnergy kineticEnergy totalEnergy" << std::endl;
    total = 0.0;

    //int velScaleEvery = 200; // while equilibrating, scale the velocities every so many turns

    int nsteps = equilSteps + prodSteps;

    ensembleAvgDiffusion = std::vector<double> (nsteps, 0.0);

    // make the file msd.dat
    std::ostringstream stringStream3;
    stringStream3.flush();
    stringStream3.str("");
    stringStream3 << name << "_msd.dat";
    std::string msdData = stringStream3.str();
    std::ofstream msdDataFile(msdData.c_str(), std::ios::out);

    msdDataFile << "# format: step LJTime meanSquareDisplacement diffusionCoeff" << std::endl;

    // get the current potential energy of the system; also initializes the tail correction
    
    ComputeTotalPotentialEnergy();

    // we are initializing the simulation; get the current force on each atom
    //zeroForces();

    //ComputeTotalForce();

    // have this call 'zeroForces' first
    calculateForces();
    // if we are doing our equilibration, we can now write the forces data to file
    // -- note that we expect the net force on each atom to be zero
    if (!(production)) {
        std::ostringstream stringStream2;
        stringStream2.flush();
        stringStream2.str("");
        stringStream2 << name << "_init_forces.dat";
        std::string initForces = stringStream2.str();
        std::ofstream initForcesFile(initForces.c_str(), std::ios::out);
        initForcesFile << "# format: atomIndex Fx Fy Fz Vx Vy Vz" << std::endl;
        initForcesFile << atoms.size() << "\n" << std::endl;
        std::vector<double> forcesToPrint = std::vector<double> (3, 0.0);
        std::vector<double> velsToPrint = std::vector<double> (3, 0.0);
        for (int ijk = 0; ijk < (int) atoms.size(); ijk++) {
            forcesToPrint = atoms[ijk].getForces();
            velsToPrint = atoms[ijk].getVelocities();
            initForcesFile << ijk << "      " << forcesToPrint[0] << "      " << forcesToPrint[1] << "      " 
                << forcesToPrint[2] << std::endl;
                
                /*
                "      " << velsToPrint[0] << "      " << velsToPrint[1]   << "      " 
                << velsToPrint[2] << std::endl;
                */
        };
    };

    //double testDistance;
    
    //std::cout << atoms[0].xs[0] << " " << atoms[0].xs[1] << " " << atoms[0].xs[2] << std::endl;
    //std::cout << atoms[215].xs[0] << " " << atoms[215].xs[1] << " " << atoms[215].xs[2] << std::endl;

    blockSumKinetic = 0.0;
    //exit(1);

    for (int i = 0; i < nsteps; i++) {
        // let's just put all the simulation into one 'run' call for now

        if (i > equilSteps) {
            production = true;
        };
        /* Velocity-Verlet Integrator
         * --update velocities by half step, according to current forces (accelerations)
         * --update positions, according to current velocities & accelerations
         * --zero the forces
         * --compute the forces present in the system
         * --update the velocities by a half step
         * ----and do this 'nsteps' number of times
         */

        // half-step velocities update
        velocityHalfStep();

        // full step positions update
        // --note: apply PBC's within this for each atom as it is translated!
        positionStep();

        //isokinetic rescaling of velocities during equilibration
        /*
        if (!(production)) {
            if (i%velScaleEvery == 0 && (i != 0)) {
                rescaleVelocities();
            }
        }
        */
        // zero the fs vector for each atom
        zeroForces();

        // re-calculate the forces
        calculateForces();

        // update the velocities by a halfstep
        velocityHalfStep();

        // take the average of the deltaR2 vector right now
        calculateAvgDeltaR2(i);

        ComputeTotalPotentialEnergy();
        
        //std::cout << "Potential Energy at step " << i << " is " << totalPE << std::endl;

        ComputeTotalKineticEnergy();
        // if we are at the end of printEvery during a production cycle, do stuff:
        blockSumKinetic += totalKE;


        if ( (i % 200) == 0) {
            std::cout << "At step " << i << std::endl;
            //printConfig(name,i);
        };
        if ( ( (i % printXYZ ) == 0 ) and (production) ) {
            //std::cout << "At step " << i << " of production." << std::endl;
            printConfig(name,i);           
        };

        // output to simData the following quantities each step:
        // step   PE    KE   (totalE)     
        double totalPE_n, totalE_n;
        // we print our properties out normalized s.t. they are per-atom
        totalPE_n = totalPE / ((double) atoms.size());
        totalKE_n = totalKE / ((double) atoms.size());
        totalE_n = totalPE_n + totalKE_n;

        simDataFile << i  << " " << totalPE_n << " " << totalKE_n << " " << (totalE_n) << std::endl;

        double time = i * dt / (sqrt( sigma * sigma * mass / epsilon));
        double diffusion;
        if (time > 0.0) {
            diffusion = ensembleAvgDiffusion[i] / (6.0 * time);
        msdDataFile << i << " " << time << " " << ensembleAvgDiffusion[i] << " " << diffusion << " " << std::endl;
        };
    };
}

void Simulation::velocityHalfStep() {
    // here, we modify the velocities by a halfstep according to the velocity-verlet integration scheme
    std::vector<double> dv = std::vector<double> (3, 0.0);;
    std::vector<double> thisForce; 
    // NOTE: at this point, each atom should know all pairwise forces acting on it;
    // they just need to be integrated s.t. the velocity is modified accordingly
    for (int i = 0; i < (int) atoms.size(); i++) {
        Atom &thisAtom = atoms[i];
        thisForce = thisAtom.getForces();
        dv[0] = dt2m * thisForce[0];
        dv[1] = dt2m * thisForce[1];
        dv[2] = dt2m * thisForce[2];

        thisAtom.addVelocities(dv);
    };
}

void Simulation::positionStep() {
    std::vector<double> vel;
    std::vector<double> dx = std::vector<double> (3,0.0);
    double dr2; // deltaR^2
    // NOTE: at this point, each atom should know its velocity;
    // we simply need to integrate it s.t. the position is modified accordingly
    // ALSO, here, we might keep track of the < ( delta r)^2 > quantity;
    // additionally, we would enforce PBC in this function /after/ modifying the position

    for (int i = 0; i < (int) atoms.size(); i++) {
        Atom &thisAtom = atoms[i];
        vel = thisAtom.getVelocities();
        dx[0] = dt * vel[0];
        dx[1] = dt * vel[1];
        dx[2] = dt * vel[2];

        dr2 = (dx[0] * dx[0]) + (dx[1] * dx[1]) + (dx[2] * dx[2]);

        deltaR2[i] += dr2;

        thisAtom.addCoordinates(dx);

        box->enforcePBC(thisAtom);
    };
}


void Simulation::rescaleVelocities() {
    double blockCount = 200.0;
    double nAtoms = (double) atoms.size();
    blockSumKinetic /= (blockCount * nAtoms);
    std::vector<double> theseVels = std::vector<double> (3, 0.0);

    for (int i = 0; i < (int) atoms.size(); i++) {
        theseVels = atoms[i].getVelocities();
        theseVels[0] *= (sqrt(3.0 * Tstar / (2.0 * blockSumKinetic)));
        theseVels[1] *= (sqrt(3.0 * Tstar / (2.0 * blockSumKinetic)));
        theseVels[2] *= (sqrt(3.0 * Tstar / (2.0 * blockSumKinetic)));
        atoms[i].setVelocities(theseVels[0],theseVels[1],theseVels[2]);
    };
    blockSumKinetic = 0.0;
}

void Simulation::calculateAvgDeltaR2(int thisStep) {
    // here we get the diffusion coefficient as a function of timestep
    double thisSum = 0.0;
    for (int i = 0; i < (int) atoms.size(); i++) {
        thisSum += box->computeAbsoluteDistanceTraveled(atoms[i]);
    };

    thisSum /= ( (double) atoms.size());
    ensembleAvgDiffusion[thisStep] = thisSum;
}

void Simulation::ComputeTotalKineticEnergy() {
    double thisSum = 0.0;
    totalKE = 0.0;
    std::vector<double> vels;
    for (int i = 0; i < (int) atoms.size(); i++) {
        vels = atoms[i].getVelocities();
        thisSum += ((1.0 / 2.0) * (mass) * ((vels[0] * vels[0]) + 
                                            (vels[1] * vels[1]) + 
                                            (vels[2] * vels[2])) );
    };
    totalKE = thisSum;
}

void Simulation::zeroForces() {
    // set the fs vector to (0.0, 0.0, 0.0) for all atoms in simulation
    for (int i = 0; i < (int) atoms.size(); i++) {
        Atom &a = atoms[i];
        a.resetForces();
    };
}

// update fs vector for each atom
void Simulation::calculateForces() {
    zeroForces();
    double thisForce = 0.0;
    //double distance = 0.0;
    for (int i = 0; i < (int) (atoms.size() - 1) ; i++) {
        Atom &a = atoms[i];
        for (int j = i+1; j < (int) atoms.size(); j++) {
            Atom &b = atoms[j];
            thisForce = LJForce(a,b);
            //std::cout << "atoms " << i << " and " << j << " dist " << distance << " force " << thisForce << std::endl;
        };
    };
  
}

void Simulation::ComputeTotalPotentialEnergy() {
    // zero our totalPE variable
    totalPE = 0.0;
    tail = 0.0;
    for (int i = 0; i < ((int) atoms.size()- 1); i++) {
        for (int j = i+1; j < (int) atoms.size(); j++) {
            totalPE += LJPotential(i,j);
        };
    };

    // and our tail correction, since we truncate the potential abruptly
    // note that 'density' here is rhoStar - so, divide by sig3
    tail = ( (8.0 / 3.0) * M_PI * density) ;
    tail *= ( ( (1.0 / 3.0) * (1.0 / rcut9)) - (1.0 / rcut3));
    totalPE += tail;
}

double Simulation::LJPotential(int ii,int jj) {
    Atom &a = atoms[ii];
    Atom &b = atoms[jj];
    std::vector<double> rij_vector = std::vector<double> (3, 0.0);
    
    rij_vector = box->computeDistanceVector(a,b);
    
    double rij = (sqrt( rij_vector[0] * rij_vector[0] +
                        rij_vector[1] * rij_vector[1] + 
                        rij_vector[2] * rij_vector[2]));
    if (rij < rcut) {
        // we have global values for sigma, eps, sig6, sig12...

        double ijEnergy = 0.0;

        double srsq = sig2 / (rij * rij);
        double sr6 = srsq * srsq * srsq;
        double sr12 = sr6 * sr6;

        ijEnergy = 4.0 * epsilon * ( sr12 - sr6 ); 

        //std::cout << "sig2 " << sig2 << " rij " << rij << " sig2 / rij * rij " << srsq << std::endl;
        //std::cout << "srsq " << srsq << " sr6 " << sr6 << " sr12 " << sr12 << " epsilon " << epsilon << std::endl;
        //std::cout << "atoms " << ii << " and " << jj << " distance " << rij << " potential " << ijEnergy << std::endl;
        return ijEnergy;
    } else {
        return 0.0;
    };

}

double Simulation::LJForce(Atom &a, Atom &b) {
    std::vector<double> rij_vector = std::vector<double> (3, 0.0);
    
    // returns un-normalized vector going between the minimum images of i & j
    rij_vector = box->computeDistanceVector(a,b);

    std::vector<double> force_vec = std::vector<double> (3, 0.0);
    std::vector<double> force_vec_neg = std::vector<double> (3, 0.0);
    // we will need the unit vector
    double dx = 0.0;
    double dy = 0.0;
    double dz = 0.0;
    
    
    double rij = (sqrt( rij_vector[0] * rij_vector[0] +
                        rij_vector[1] * rij_vector[1] + 
                        rij_vector[2] * rij_vector[2]));

    if (rij < rcut) {
        
        double sr2 = sig2 /  (rij * rij);
        double sr6 = sr2 * sr2 * sr2;
        double sr12 = sr6 * sr6;
        
        // we will calculate the magnitude of the force and dot it with rij_unit_vec
        // to get its directionality

        // magnitude of our force from atom i to atom j
        double ijForce = (24.0 * epsilon / (rij*rij)) * ( (2.0 * sr12) - sr6 ) ;

        //std::cout << "ijForce calculated to be : " << ijForce << std::endl;

        dx = rij_vector[0] / rij;
        dy = rij_vector[1] / rij;
        dz = rij_vector[2] / rij;

        force_vec[0] = ijForce * dx;
        force_vec[1] = ijForce * dy;
        force_vec[2] = ijForce * dz;

        force_vec_neg[0] = 0.0 - force_vec[0];
        force_vec_neg[1] = 0.0 - force_vec[1];
        force_vec_neg[2] = 0.0 - force_vec[2];

        a.addForce(force_vec[0], force_vec[1], force_vec[2]);
        b.addForce(force_vec_neg[0], force_vec_neg[1], force_vec_neg[2]);

        return ijForce;

    } else {
        b.addForce(force_vec[0], force_vec[1], force_vec[2]);
        a.addForce(force_vec[0], force_vec[1], force_vec[2]);

        return 0.0;
    
    };
    
    return 0.0;

}

void Simulation::printConfig(std::string name, int step) {
    std::ostringstream stringStream;
    stringStream.flush();
    stringStream.str("");
    stringStream << name << "_step_" << step << ".xyz";
    std::string xyzFileName = stringStream.str();
    std::ofstream xyzFile(xyzFileName.c_str(), std::ios::out);
    
    // write the number of atoms, and then a comment line denoting the format

    xyzFile << atoms.size() << std::endl;

    if (step == 0) {
        xyzFile << "# format: xpos ypos zpos xvel yvel zvel" << std::endl;
    } else { 
        xyzFile << std::endl;
    };

    // now, write atom name, xpos, ypos, zpos (newline) to file
    std::vector<double> theseCoords = std::vector<double> (3,0.0);
    std::vector<double> theseVels = std::vector<double> (3, 0.0);
    for (unsigned int i = 0; i < atoms.size(); i++) {
        theseCoords = atoms[i].getCoordinates();
        theseVels = atoms[i].getVelocities();
        xyzFile << "atom" << i << " " << theseCoords[0] << " " << theseCoords[1] << " " << theseCoords[2] << std::endl;
           // << " " << theseVels[0] << " " << theseVels[1] << " " << theseVels[2] << std::endl;
    };

}


