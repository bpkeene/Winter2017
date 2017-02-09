#include "atom.h"
#include <vector>
#include <iostream>

// atom constructor
Atom::Atom() {
    // initialize the coords as a vector of length 3
    coords = std::vector<double> (3, 0.0);

    oldCoords = std::vector<double> (3, 0.0);
    //force = std::vector<double> (3, 0.0);
    ijForce = std::vector<double> (3, 0.0);

    totalForce = std::vector<double> (3, 0.0);
};

void Atom::setLJ(double _sigma, double _epsilon) {
    sigma = _sigma;
    epsilon = _epsilon;

    sig2 = sigma * sigma;
    sig6 = sig2 * sig2 * sig2;
    sig12 = sig6 * sig6;
};

void Atom::setCoordinates(double _x, double _y, double _z) {
    coords[0] = _x;
    coords[1] = _y;
    coords[2] = _z;
};

void Atom::initializeForces(int numAtoms) {
    
    forces = std::vector< std::vector<double> > (numAtoms, std::vector<double> (3,0.0));
    oldForces = std::vector< std::vector<double> > (numAtoms, std::vector<double> (3, 0.0));
};

void Atom::saveForce(int idx) {
    // atom idx is being moved; save its old force, in case the move is rejected
    ijForce[0] = forces[idx][0];
    ijForce[1] = forces[idx][1];
    ijForce[2] = forces[idx][2];
    idxSaved = idx;
};

std::vector< std::vector<double> > Atom::getOldForces() {
    return oldForces;
};

std::vector< std::vector<double> > Atom::getForcesMatrix() {
    return forces;
};

void Atom::resetForce() {
    forces[idxSaved][0] = ijForce[0];
    forces[idxSaved][1] = ijForce[1];
    forces[idxSaved][2] = ijForce[2];
};

void Atom::saveForces() {
    oldForces = forces;
};

void Atom::resetForces() {
    forces = oldForces;
};

void Atom::setForce(int idx, std::vector<double> newForce) {
    // calling this automatically saves the forces
    //std::cout << "About to call saveForce(idx) " << std::endl;
    saveForce(idx);

    //std::cout << "Setting forces[idx] = newForce" << std::endl;
    forces[idx] = newForce;
};

void Atom::computeTotalForce() {
    // sum across the forces vector to get Fx, Fy, Fz for this atom
    totalForce = std::vector<double> (3, 0.0);

    for (int i = 0; i < (int) forces[0].size() ; i++) {
        for (int j = 0; j < 3; j++ ) {
            totalForce[j] += forces[i][j];
        };
    };
};

std::vector<double> Atom::getTotalForce() {
    return totalForce;
};

std::vector<double> Atom::getForces(int idx) {
    std::vector<double> forcesToReturn = std::vector<double> (3, 0.0);
    forcesToReturn[0] = forces[idx][0];
    forcesToReturn[1] = forces[idx][1];
    forcesToReturn[2] = forces[idx][2];
    return forcesToReturn;
};

void Atom::updateTotalForce(int idx) {
    // the move was accepted; for each atom in atoms, update totalForce counter w.r.t this atom
    totalForce[0] += forces[idx][0] - ijForce[0];
    totalForce[1] += forces[idx][1] - ijForce[1];
    totalForce[2] += forces[idx][2] - ijForce[2];
};


void Atom::saveCoordinates() {
    oldCoords = coords;
};

void Atom::resetCoordinates() {
    coords = oldCoords;
};

std::vector<double> Atom::getOldCoordinates() {
    return oldCoords;
};

double Atom::getSigma() {
    return sigma;
};
double Atom::getSig2() {
    return sig2;
};
double Atom::getSig6() {
    return sig6;
};
double Atom::getSig12() {
    return sig12;
};
double Atom::getEpsilon() {
    return epsilon;
};
std::vector<double> Atom::getCoordinates() {
    return coords;
};
