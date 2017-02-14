#ifndef ATOM_H
#define ATOM_H

#include <vector>

class Atom {
    private:
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
        
    public:
        // just use a default constructor
        Atom();
        
        // set the LJ parameters
        void setLJ(double _sigma, double _epsilon);
       
        // retrieve the LJ parameters
        double getSigma();
        double getEpsilon();

        // set the coordinates (for initializations)
        void setCoordinates(double _x, double _y, double _z);
        std::vector<double> getCoordinates();
        void addCoordinates(std::vector<double> &_xs);
     
        // populates xs_init vector
        void setInitCoordinates();
        std::vector<double> getInitCoordinates();

        // during the middle of an integration step, we want to reset the forces;

        // add the forces that were calculated to the tally
        void resetForces();
        std::vector<double> getForces();
        void addForce(double &, double &, double &);

        // set the velocities (for initializations)
        void setVelocities(double _vx, double _vy, double _vz);
        std::vector<double> getVelocities();
        void addVelocities(std::vector<double> &_vs);


        void setMass(double _mass);

        std::vector<double> getBoxesTraveled();
        void addBoxesTraveled(double, double, double);
};



#endif /* ATOM_H */
