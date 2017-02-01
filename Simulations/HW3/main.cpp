//#include "simulation.h" // probably to run a simulation
#include <iostream> // for the std::cout
#include <sstream> // for std::stringstream
#include <string>
#include <cstdio>
#include <fstream>

int main(int argc, char *argv[]) {

    std::string filename;
    
    if (argc != 2) {
        std::cout << "No input file provided, or too many arguments.  Exiting." << std::endl;
        return(1);
    }; 
    
    if ((argv[1] == NULL)) {
        return(1);
    };

    std::ifstream fileStream(argv[1]);
    
    if (fileStream.is_open()) {
        char x;
        while (fileStream.get(x)) {
               std::cout << x;
        };
    };

};













