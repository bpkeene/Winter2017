#ifndef ATOM_H
#define ATOM_H

#include <vector>

class Atom {
    public:

        // every LJ atom will have a sigma, epsilon, and tuple of coordinates;
        double sigma;
        double epsilon;
        double mass;  // need mass for MD simulation
        
        //position, force, and velocity vectors for this atom
        std::vector<double> xs;
        std::vector<double> fs;
        std::vector<double> vs;      

        // a vector of our initial coords
        std::vector<double> xs_init;

        // a vector of our initial velocities
        std::vector<double> vs_init;

        // a vector of 3 int's, initialized to (0,0,0), denoting lengths 
        // in the x, y, z direction traveled since  beginning of the simulation
        std::vector<double> boxesTraveled;       
        
        // just use a default constructor
        Atom();
        
        // set the LJ parameters
        void setLJ(double _sigma, double _epsilon);
       
        // retrieve the LJ parameters
        double getSigma();
        double getEpsilon();

        // set the coordinates (for initializations)
        void setCoordinates(double _x, double _y, double _z);
     
        // populates xs_init vector
        void setInitCoordinates();

        // during the middle of an integration step, we want to reset the forces;
        void resetForces();

        // add the forces that were calculated to the tally
        void addForce(std::vector<double> _fs);

        // set the velocities (for initializations)
        void setVelocities(double _vx, double _vy, double _vz);

        void addVelocities(std::vector<double> _vs);

        void addCoordinates(std::vector<double> _xs);

        // get the coordinates
        std::vector<double> getCoordinates();

        // get the velocities
        std::vector<double> getVelocities();

        // get the forces
        std::vector<double> getForces();

        void setMass(double _mass);
};



#endif /* ATOM_H */
