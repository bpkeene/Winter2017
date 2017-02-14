#ifndef BOX_H
#define BOX_H

#include "atom.h"
#include "random_mars.h"
#include <vector>


class Box {
    private:
        double x_dim;
        double y_dim;
        double z_dim;

        double sigma;
        double volume;

        std::vector<double> dimensions;
        
    public:
        Box();
        // our constructor
        Box(int nAtoms, double rhoStar, double sigma);
        
        void initializeAtoms(std::vector<Atom> &, double);

        // get the minimum image distance between two atoms
        double computeDistance(Atom &, Atom &);

        // calculate the absolute distance ( (delta r) ^ 2) traveled from the initial position
        double computeAbsoluteDistanceTraveled(Atom &);

        std::vector<double> getDim();

        double getVolume();

        void enforcePBC(Atom &);
};

#endif /* BOX_H */
