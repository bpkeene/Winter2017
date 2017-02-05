#ifndef BOX_H
#define BOX_H

#include "atom.h"
#include <vector>


class Box {
    private:
        double x_dim;
        double y_dim;
        double z_dim;

        double volume;

        std::vector<double> dimensions;
        
    public:
        Box();
        // our constructor
        Box(int nAtoms, double rhoStar, double sigma);

        double computeDistance(Atom &, Atom &);
        std::vector<double> getDim();
};

#endif /* BOX_H */
