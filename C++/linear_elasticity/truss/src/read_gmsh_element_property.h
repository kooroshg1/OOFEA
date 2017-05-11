//
// Created by koorosh on 12/7/16.
//
#include <iostream>
#include <map>
#include <vector>
#ifndef TEST15_READ_TRUSS_ELEMENT_PROPERTIES
#define TEST15_READ_TRUSS_ELEMENT_PROPERTIES
namespace Extra{
    class elementProperty {
    public:
        elementProperty();
        void read(std::string file_name);
        
        ~elementProperty();
        
        std::map<unsigned int, std::vector<double>> material_property;
        std::map<unsigned int, std::vector<double>> cross_section_property;
    };
}
#endif //TEST15_READ_TRUSS_ELEMENT_PROPERTIES
