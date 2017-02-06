#include "atom.h"
#include <vector>

// atom constructor
Atom::Atom() {
    // initialize the coords as a vector of length 3
    coords = std::vector<double> (3, 0.0);

    oldCoords = std::vector<double> (3, 0.0);
    // initialize the forces as a vector of length 3
    force = std::vector<double> (3, 0.0);

};

void Atom::setLJ(double _sigma, double _epsilon) {
    sigma = _sigma;
    epsilon = _epsilon;

};

void Atom::setCoordinates(double _x, double _y, double _z) {
    coords[0] = _x;
    coords[1] = _y;
    coords[2] = _z;
};

void Atom::saveCoordinates() {
    oldCoords = coords;
};

void Atom::resetCoordinates() {
    coords = oldCoords;
};

double Atom::getSigma() {
    return sigma;
};
double Atom::getEpsilon() {
    return epsilon;
};
std::vector<double> Atom::getCoordinates() {
    return coords;
};
