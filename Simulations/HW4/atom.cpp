#include "atom.h"
#include <vector>
#include <iostream>

// atom constructor
Atom::Atom() {
    // initialize vectors of length 3
    xs = std::vector<double> (3, 0.0);

    fs = std::vector<double> (3, 0.0);

    vs = std::vector<double> (3, 0.0);

    boxesTraveled = std::vector<double> (3, 0.0);

    xs_init = std::vector<double> (3, 0.0);
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
std::vector<double> Atom::getInitCoordinates() {
    return xs_init;
};

void Atom::addCoordinates(std::vector<double> &_xs) {
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
    fs = std::vector<double> (3, 0.0);
};

void Atom::addForce(double &fsx, double &fsy, double &fsz) {
    fs[0] += fsx;
    fs[1] += fsy;
    fs[2] += fsz;
};

std::vector<double> Atom::getForces() {
    return fs;
};
void Atom::addVelocities(std::vector<double> &_vs) {
    vs[0] += _vs[0];
    vs[1] += _vs[1];
    vs[2] += _vs[2];
};

void Atom::setVelocities(double _vx,double _vy,double _vz) {
    vs[0] = _vx;
    vs[1] = _vy;
    vs[2] = _vz;
};

std::vector<double> Atom::getVelocities() {
    return vs;
};

std::vector<double> Atom::getBoxesTraveled() {
    return boxesTraveled;
};

void Atom::addBoxesTraveled(double _x, double _y, double _z) {
    boxesTraveled[0] += _x;
    boxesTraveled[1] += _y;
    boxesTraveled[2] += _z;
};

std::vector<double> Atom::getCoordinates() {
    return xs;
};
