#include "atom.h"
#include <vector>
#include <iostream>

// atom constructor
Atom::Atom() {
    // initialize vectors of length 3
    xs = std::vector<double> (3, 0.0);

    fs = std::vector<double> (3, 0.0);

    vs = std::vector<double> (3, 0.0);

    boxesTraveled = std::vector<double> (3, 0);
};

void Atom::setMass(double _mass) {
    mass = _mass;
};

void Atom::setLJ(double _sigma, double _epsilon) {
    sigma = _sigma;
    epsilon = _epsilon;
};

void Atom::setCoordinates(double _x, double _y, double _z) {
    xs[0] = _x;
    xs[1] = _y;
    xs[2] = _z;
};

// this should /only/ be called when box initializes the atoms!
void Atom::setInitCoordinates() {
    xs_init = xs;
};

void Atom::addCoordinates(std::vector<double> _xs) {
    xs[0] += _xs[0];
    xs[1] += _xs[1];
    xs[2] += _xs[2];
};

double Atom::getSigma() {
    return sigma;
};
double Atom::getEpsilon() {
    return epsilon;
};

void Atom::resetForces() {
    fs = std::vector<double (3, 0.0);
};

void Atom::addForce(std::vector<double> _fs) {
    fs[0] += _fs[0];
    fs[1] += _fs[1];
    fs[2] += _fs[2];
};

void Atom::addVelocities(std::vector<double> _vs) {
    vs[0] += _vs[0];
    vs[1] += _vs[1];
    vs[2] += _vs[2];
};

void Atom::setVelocities(std::vector<double> _vs) {
    vs[0] = _vs[0];
    vs[1] = _vs[1];
    vs[2] = _vs[2];
};


std::vector<double> Atom::getCoordinates() {
    return xs;
};
