#ifndef ATOM_H
#define ATOM_H

#include <vector>

class Atom {
    private:

        // every LJ atom will have a sigma, epsilon, and tuple of coordinates;
        double sigma;
        double epsilon;
        std::vector<double> coords;
        std::vector<double> force;
        std::vector<double> oldCoords;
    public:
        Atom();
        
        // set the LJ parameters
        void setLJ(double _sigma, double _epsilon);
       
        // retrieve the LJ parameters
        double getSigma();
        double getEpsilon();

        // set the coordinates
        void setCoordinates(double _x, double _y, double _z);

        void saveCoordinates();

        void resetCoordinates();

        // get the coordinates
        std::vector<double> getCoordinates();
};



#endif /* ATOM_H */
