//
// Created by koorosh on 12/7/16.
//
#include "read_gmsh_bc.h"
#include <vector>
#include <map>
#include <iterator>
#include <fstream>
#include <iostream>

namespace Extra {
    boundaryCondition::boundaryCondition() {

    }

    boundaryCondition::~boundaryCondition() {

    }

    void boundaryCondition::read(std::string file_name) {
        std::string file_line; // For reading file line by line
        std::ifstream file_address(file_name.c_str());

        unsigned int number_of_dirichlet_bc;
        unsigned int number_of_force_bc;
        while (true) {
            std::getline(file_address, file_line);
            // Add displacement boundary to a std::map container
            if (file_line.find("$DisplacementBoundary") == static_cast<std::string::size_type>(0)) {
                file_address >> number_of_dirichlet_bc;
                for (unsigned int i=0; i != number_of_dirichlet_bc; i++) {
                    unsigned int node_number, degree_of_freedom;
                    file_address >> node_number >> degree_of_freedom;
                    // Since libMesh numbering starts at 0,
                    // I need to subtract the node_number by 1.
                    displacement_boundary.insert(std::make_pair(node_number - 1, degree_of_freedom));
                }
            }

            // Add force boundary to a std::map container
            if (file_line.find("$ForceBoundary") == static_cast<std::string::size_type>(0)) {
                file_address >> number_of_force_bc;
                for (unsigned int i=0; i != number_of_force_bc; i++) {
                    unsigned int node_number;
                    double f_x, f_y, f_z;
                    std::vector<double> load;
                    file_address >> node_number >> f_x >> f_y >> f_z;
                    load.push_back(f_x);
                    load.push_back(f_y);
                    load.push_back(f_z);
                    // I use this to check if the load is NOT applied
                    // at the corresponding node. '1' means load is
                    // NOT applied and '0' means load IS applied.
                    load.push_back(1);
                    // Since libMesh numbering starts at 0,
                    // I need to subtract the node_number by 1.
                    load_boundary.insert(std::make_pair(node_number - 1, load));
                }
            }

            if (file_address.eof()) {
                break;
            }
        }
//        std::map<unsigned int, unsigned int>::iterator it = displacement_boundary.begin();
//        while (it != displacement_boundary.end()) {
//            std::cout << it->first << "\t" << it->second << std::endl;
//            it++;
//        }
//
//        std::map<unsigned int, std::vector<double>>::iterator it_foce = load_boundary.begin();
//        while (it_foce != load_boundary.end()) {
//            std::cout << it_foce->first << "\t"
//                      << it_foce->second[0] << "\t"
//                      << it_foce->second[1] << "\t"
//                      << it_foce->second[2] << "\t" << std::endl;
//            it_foce++;
//        }
    }
}
