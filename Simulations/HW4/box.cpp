#include "box.h"
#include "atom.h"
#include <vector>
#include <math.h>
#include <iostream>

// our constructor
Box::Box(int nAtoms, double rhoStar, double _sigma) {
    
    sigma = _sigma;
    dimensions = std::vector<double> (3,0.0);
    volume = ( (double) nAtoms) / rhoStar;
    std::cout << "volume calculated to be: " << volume << std::endl;

    x_dim = pow(volume, (1.0/3.0));
    y_dim = x_dim;
    z_dim = x_dim;
    
    std::cout << "x, y z dims: " << x_dim << std::endl;
    dimensions[0] = x_dim;
    dimensions[1] = y_dim;
    dimensions[2] = z_dim;
}

void Box::initializeAtoms(std::vector<Atom> &atoms, double _mass) {

    double increment;
    //double initial_x_pos, initial_y_pos, initial_z_pos;
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
    //initial_x_pos = 0.0;
    //initial_y_pos = 0.0;
    //initial_z_pos = 0.0;

    xpos = ypos = zpos = 0.0;

    int atomIndex = 0;
    int numberOfAtoms = (int) atoms.size();
    double boxLen = x_dim;
    double numPerDim = (double) atomsPerLength;
    std::cout << "about to initialize the atoms!" << std::endl;
    for (int ix = 0; ix < atomsPerLength; ix++) {
        //zpos = initial_z_pos + i * increment;
        for (int iy = 0; iy < atomsPerLength; iy++) {
            //ypos = initial_y_pos + j * increment;
            //xpos = initial_x_pos;
            for (int iz = 0; iz < atomsPerLength; iz++) {
                if (atomIndex < numberOfAtoms) {
                    double x = (ix + 0.5) * boxLen / numPerDim;
                    double y = (iy + 0.5) * boxLen / numPerDim;
                    double z = (iz + 0.5) * boxLen / numPerDim;
                    atoms[atomIndex].setCoordinates(x,y,z);
                    atoms[atomIndex].setInitCoordinates();
                    atoms[atomIndex].setMass(_mass);
                    xpos += increment;
                    atomIndex++;
                } else {
                    break;
                };
            };
        };
    };
}

// for some pair of atoms, compute the distance between their minimum images
double Box::computeDistance(Atom &atom1, Atom &atom2) {
 
    double dx, dy, dz;
    double dxm = 0.0;
    double dym = 0.0;
    double dzm = 0.0;
    double ret;
    double boxlenhalf = 0.5 * x_dim;
    // get the coords of the atoms
    std::vector<double> coords1 = atom1.getCoordinates();
    std::vector<double> coords2 = atom2.getCoordinates();

    // get dx, dy, dz
    dx = coords1[0] - coords2[0];
    dy = coords1[1] - coords2[1];
    dz = coords1[2] - coords2[2];

    if (dx > boxlenhalf) dxm = dxm - x_dim;
    if (dx < -boxlenhalf) dxm = dxm + x_dim;
    if (dy > boxlenhalf) dym = dym - y_dim;
    if (dy < -boxlenhalf) dym = dym + y_dim;
    if (dz > boxlenhalf) dzm = dzm - z_dim;
    if (dz < -boxlenhalf) dzm = dzm + z_dim;

    dx += dxm;
    dy += dym;
    dz += dzm;

    // TODO: should this be divided by sigma?
    ret = (sqrt((dx*dx) + (dy*dy) + (dz*dz)));

    /*
    // get the minimum image dx, dy, dz
    dxm = fabs(fabs(dx) - x_dim);
    dxm = std::min(fabs(dx), dxm);

    dym = fabs(fabs(dy) - y_dim);
    dym = std::min(fabs(dy), dym);

    dzm = fabs(fabs(dz) - z_dim);
    dzm = std::min(fabs(dz), dzm);

    */
    // return the distance
    return ret;

}

void Box::enforcePBC(Atom &atom) {
    std::vector<double> theseCoords = atom.getCoordinates();
    std::vector<double> boxesTraveled = std::vector<double> (3, 0.0);
    // apply PBC
    if (theseCoords[0] < 0.0) {
        theseCoords[0] += x_dim;
        boxesTraveled[0] -= 1.0;
    };
    if (theseCoords[0] > x_dim) {
        theseCoords[0] -= x_dim;
        boxesTraveled[0] += 1.0;
    };
    if (theseCoords[1] <0.0) {
        theseCoords[1] += y_dim;
        boxesTraveled[1] -= 1.0;
    };
    if (theseCoords[1] > y_dim) {
        theseCoords[1] -= y_dim;
        boxesTraveled[1] += 1.0;
    };
    if (theseCoords[2] < 0.0) {
        theseCoords[2] += z_dim;
        boxesTraveled[2] -= 1.0;
    };
    if (theseCoords[2] > z_dim) {
        theseCoords[2] -= z_dim;
        boxesTraveled[2] += 1.0;
    };

    // could probably make it so that setCoordinates takes a vector.. but eh
    atom.setCoordinates(theseCoords[0], theseCoords[1], theseCoords[2]);
    atom.addBoxesTraveled(boxesTraveled[0], boxesTraveled[1], boxesTraveled[2]);
}

double Box::computeAbsoluteDistanceTraveled(Atom &atom) {
    // compute the absolute distance a particle has traversed, unwrapping the PBC
    std::vector<double> coords = atom.getCoordinates();

    // copy the initial values, so we don't overwrite them
    std::vector<double> coords_init = atom.getInitCoordinates();

    // and finally, copy the boxes traveled vector;
    std::vector<double> boxes = atom.getBoxesTraveled();

    double deltar2;
    double dx, dy, dz;
    
    // the second term (e.g., "x_dim * boxes[0]") amounts to unwrapping the periodicity
    // which permits us to get the absolute quantity that we desire
    dx = (coords[0] + (x_dim * boxes[0])) - coords_init[0];
    dy = (coords[1] + (y_dim * boxes[1])) - coords_init[1];
    dz = (coords[2] + (z_dim * boxes[2])) - coords_init[2];

    deltar2 = (dx * dx) + (dy * dy) + (dz * dz);
    return deltar2;

}
double Box::getVolume() {
    return volume;
}


// return the lengths of the sides of the box
std::vector<double> Box::getDim() {
    return dimensions;
}


