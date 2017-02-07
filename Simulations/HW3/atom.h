#ifndef ATOM_H
#define ATOM_H

#include <vector>

class Atom {
    private:

        // every LJ atom will have a sigma, epsilon, and tuple of coordinates;
        double sigma;
        double epsilon;
        double sig2;
        double sig6;
        double sig12;
        std::vector<double> coords;
        std::vector<double> oldCoords;
        std::vector< std::vector<double> > forces;

        // a 3-value vector containing the sum of all forces on this atom 
        // at this point in time;
        // calling this computes the force using the forces vector
        // and so the forces vector /must/ be up to date!
        std::vector<double> totalForce;

        // this atom is being moved; store all the old forces
        std::vector< std::vector<double> > oldForces;
       
        // another atom is being moved; store the old ij force
        std::vector<double> ijForce;
        int idxSaved; // index to which ijForce corresponds

    public:
        Atom();
        
        // set the LJ parameters
        void setLJ(double _sigma, double _epsilon);
       
        // retrieve the LJ parameters
        double getSigma();

        double getEpsilon();

        double getSig2();
        double getSig6();
        double getSig12();

        void initializeForces(int); // initializes the force vectors, nothing more;
                // does /not/ populate them!
        
        void saveForce(int); // save forces of atom w.r.t. atom idx
        void saveForces(); // save /all/ forces; this atom is being moved!
        std::vector< std::vector<double> > getOldForces(); // get the old forces that were saved
        std::vector< std::vector<double> > getForcesMatrix();
        std::vector<double> getOldCoordinates();
        // before this atom was perturbed
        void setForce(int idx, std::vector<double>);
        
        void resetForce(); // reset forces w.r.t. atom idxSaved
        void resetForces(); // reset /all/ forces (this atom failed to move!)
        
        void computeTotalForce(); // populate totalForce vector
        void updateTotalForce(int); // update the totalForce vector after each move
        std::vector<double> getTotalForce(); // retrieve the total force vector

        // set the coordinates
        void setCoordinates(double _x, double _y, double _z);
        
        void saveCoordinates();

        void resetCoordinates();

        // get some subset of the vector forces above;
        std::vector<double>  getForces(int);

        // get the coordinates
        std::vector<double> getCoordinates();
};



#endif /* ATOM_H */
