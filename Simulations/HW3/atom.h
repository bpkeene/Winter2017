#ifndef ATOM_H
#define ATOM_H


class Atom {
    private:

        // every LJ atom will have a sigma, epsilon, and tuple of coordinates;
        double sigma;
        double epsilon;
        std::vector<double> coords;

    public:
        Atom() {};
        
        // set the LJ parameters
        void setLJ(double, double);
        
        // set the coordinates
        void setCoordinates(double, double, double);

        // get the coordinates
        std::vector<double> getCoordinates();
};



#endif /* ATOM_H */
