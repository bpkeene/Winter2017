#include "box.h"
#include "atom.h"
#include <vector>
#include <math.h>

// our constructor
Box::Box(int nAtoms, double rhoStar, double sigma) {
    dimensions = std::vector<double> (3,0.0);
    volume = (sigma * sigma * sigma) * nAtoms / rhoStar;
    
    x_dim = pow(volume, (1.0/3.0));
    y_dim = x_dim;
    z_dim = x_dim;

    dimensions[0] = x_dim;
    dimensions[1] = y_dim;
    dimensions[2] = z_dim;
}

// for some pair of atoms, compute the distance between their minimum images
double Box::computeDistance(Atom &atom1, Atom &atom2) {
 
    double dx, dy, dz;
    double dxm, dym, dzm;

    // get the coords of the atoms
    std::vector<double> coords1 = atom1.getCoordinates();
    std::vector<double> coords2 = atom2.getCoordinates();

    // get dx, dy, dz
    dx = coords1[0] - coords2[0];
    dy = coords1[1] - coords2[1];
    dz = coords1[2] - coords2[2];

    // get the minimum image dx, dy, dz
    dxm = fabs(fabs(dx) - x_dim);
    dxm = std::min(fabs(dx), dxm);

    dym = fabs(fabs(dy) - y_dim);
    dym = std::min(fabs(dy), dym);

    dzm = fabs(fabs(dz) - z_dim);
    dzm = std::min(fabs(dz), dzm);

    // return the distance
    return sqrt(dxm*dxm + dym*dym + dzm*dzm);

};

// return the lengths of the sides of the box
std::vector<double> Box::getDim() {
    return dimensions;
};


