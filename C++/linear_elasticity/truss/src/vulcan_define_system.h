//
// Created by Koorosh Gobal on 4/3/17.
//

#ifndef VULCAN_DEVELOPMENT_VULCAN_DEFINE_SYSTEM_H
#define VULCAN_DEVELOPMENT_VULCAN_DEFINE_SYSTEM_H

#include <string>

#include <libmesh/mesh.h>
#include <libmesh/equation_systems.h>

namespace vulcan {
    class FEA {
    public:
        FEA(libMesh::Mesh mesh, std::string type);
        void hello();
    };
}

#endif //VULCAN_DEVELOPMENT_VULCAN_DEFINE_SYSTEM_H
