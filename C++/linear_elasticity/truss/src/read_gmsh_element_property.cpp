//
// Created by koorosh on 12/7/16.
//
#include "read_gmsh_element_property.h"
#include <vector>
#include <map>
#include <iterator>
#include <fstream>
#include <iostream>

namespace Extra {
    elementProperty::elementProperty() {
        
    }
    
    elementProperty::~elementProperty() {
        
    }
    
    void elementProperty::read(std::string file_name) {
        std::string file_line; // For reading file line by line
        std::ifstream file_address(file_name.c_str());
        
        unsigned int number_of_material_properties;
        unsigned int number_of_element_properties;
        while (true) {
            std::getline(file_address, file_line);
            // Add displacement boundary to a std::map container
            if (file_line.find("$isotropicMaterial") == static_cast<std::string::size_type>(0)) {
                file_address >> number_of_material_properties;
                for (unsigned int i=0; i != number_of_material_properties; i++) {
                    unsigned int subdomain_id;
                    double elastic_modulus;
                    double density;
                    unsigned int calculate_sensitivity;
                    std::vector<double> material;
                    file_address >> subdomain_id >> elastic_modulus >> density >> calculate_sensitivity;
                    material.push_back(elastic_modulus);
                    material.push_back(density);
                    material.push_back(calculate_sensitivity);
                    material_property.insert(std::make_pair(subdomain_id, material));
                }
            }
            
            // Add force boundary to a std::map container
            if (file_line.find("$trussElement") == static_cast<std::string::size_type>(0)) {
                file_address >> number_of_element_properties;
                for (unsigned int i=0; i != number_of_element_properties; i++) {
                    unsigned int subdomain_id;
                    float elem_property;
                    unsigned int calculate_sensitivity;
                    std::vector<double> element;
                    file_address >> subdomain_id >> elem_property >> calculate_sensitivity;
                    element.push_back(elem_property);
                    element.push_back(calculate_sensitivity);
                    cross_section_property.insert(std::make_pair(subdomain_id, element));
                }
            }
            
            if (file_address.eof()) {
                break;
            }
        }
//    std::map<unsigned int, unsigned int>::iterator it = displacement_boundary.begin();
//    while (it != displacement_boundary.end()) {
//        std::cout << it->first << "\t" << it->second << std::endl;
//        it++;
//    }
//
//    std::map<unsigned int, std::vector<double>>::iterator it_foce = load_boundary.begin();
//    while (it_foce != load_boundary.end()) {
//        std::cout << it_foce->first << "\t"
//                  << it_foce->second[0] << "\t"
//                  << it_foce->second[1] << "\t"
//                  << it_foce->second[2] << "\t" << std::endl;
//        it_foce++;
//    }
    }
}