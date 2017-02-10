#include "box.h"
#include "atom.h"
#include <vector>
#include <math.h>
#include <iostream>

// our constructor
Box::Box(int nAtoms, double rhoStar, double sigma) {
    
    dimensions = std::vector<double> (3,0.0);
    volume = (sigma * sigma * sigma) * nAtoms / rhoStar;
    std::cout << "volume calculated to be: " << volume << std::endl;

    x_dim = pow(volume, (1.0/3.0));
    y_dim = x_dim;
    z_dim = x_dim;
    
    std::cout << "x, y z dims: " << x_dim << std::endl;
    dimensions[0] = x_dim;
    dimensions[1] = y_dim;
    dimensions[2] = z_dim;
}

void Box::initializeAtoms(std::vector<Atom> &atoms) {

    double increment;
    double initial_x_pos, initial_y_pos, initial_z_pos;
    double xpos, ypos, zpos;

    // cast number of atoms as double
    double nAtoms = (double) (atoms.size());

    // atoms per length is N^(1/3); put as ceiling, to make 
    // sure we can fit all the atoms
    double atomsPerLengthD = ceil(pow(nAtoms, (1.0/3.0))); 

    // get the dimensions of our box
    std::vector<double> maxDimensions = dimensions;


    // divide x-dim by atoms per length to get the increment
    // since it is cubic, this increment will be the same in all directions
    increment = maxDimensions[0] / atomsPerLengthD;
    int atomsPerLength = lrint(atomsPerLengthD);

    // we evenly distribute the atoms on a cubic lattice;
    initial_x_pos = 0.0;
    initial_y_pos = 0.0;
    initial_z_pos = 0.0;

    xpos = ypos = zpos = 0.0;

    int atomIndex = 0;
    int numberOfAtoms = (int) atoms.size();
    for (int i = 0; i < atomsPerLength; i++) {
        zpos = initial_z_pos + i * increment;
        for (int j = 0; j < atomsPerLength; j++) {
            ypos = initial_y_pos + j * increment;
            xpos = initial_x_pos;
            for (int k = 0; k < atomsPerLength; k++) {
                if (atomIndex < numberOfAtoms) {
                    atoms[atomIndex].setCoordinates(xpos,ypos,zpos);
                    xpos += increment;
                    atomIndex++;
                } else {
                    break;
                };
            };
        };
    };
};

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

void Box::enforcePBC(Atom &atom) {
    std::vector<double> coords = atom.getCoordinates();

    // apply PBC
    if (coords[0] < 0.0) {
        coords[0] += x_dim;
    };
    if (coords[0] > x_dim) {
        coords[0] -= x_dim;
    };
    if (coords[1] <0.0) {
        coords[1] += y_dim;
    };
    if (coords[1] > y_dim) {
        coords[1] -= y_dim;
    };
    if (coords[2] < 0.0) {
        coords[2] += z_dim;
    };
    if (coords[2] > z_dim) {
        coords[2] -= z_dim;
    };

    // could probably make it so that setCoordinates takes a vector.. but eh
    atom.setCoordinates(coords[0], coords[1], coords[2]);

};

double Box::getVolume() {
    return volume;
};


// return the lengths of the sides of the box
std::vector<double> Box::getDim() {
    return dimensions;
};


