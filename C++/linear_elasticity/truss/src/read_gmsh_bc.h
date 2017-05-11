//
// Created by koorosh on 12/7/16.
//
#include <iostream>
#include <map>
#include <vector>

#ifndef TEST13_READ_BC_READ_GMSH_BC_H
#define TEST13_READ_BC_READ_GMSH_BC_H
namespace Extra{
    class boundaryCondition {
    public:
        boundaryCondition();
        void read(std::string file_name);
        
        ~boundaryCondition();
        
        std::map<unsigned int, unsigned int> displacement_boundary;
        std::map<unsigned int, std::vector<double>> load_boundary;
    };
}
#endif //TEST13_READ_BC_READ_GMSH_BC_H
