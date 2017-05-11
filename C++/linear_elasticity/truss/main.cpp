// System includes
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <map>
#include <iterator>
#include <sstream>

// libMesh includes
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/gmsh_io.h"
#include "libmesh/equation_systems.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/dof_map.h"
#include "libmesh/elem.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/dense_submatrix.h"
#include "libmesh/node.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/vtk_io.h"

// Extra includes
#include "src/read_gmsh_element_property.h"
#include "src/read_gmsh_bc.h"

// Function prototypes
void add_structural_variables(libMesh::EquationSystems & equation_system,
                              const std::string & system_name);
void assemble_system(libMesh::EquationSystems & equation_system,
                     const std::string & system_name,
                     std::string element_properties_path);
void calculateElementOrientation(const libMesh::Elem * element,
                                 libMesh::DenseMatrix<libMesh::Real> & Re);
//void update_load(libMesh::EquationSystems & equation_system,
//                 const std::string & system_name);
void calculateElementMatrices(std::string element_type,
                              libMesh::DenseMatrix<libMesh::Number> &Ke, libMesh::DenseMatrix<libMesh::Number> &Ce, libMesh::DenseMatrix<libMesh::Number>&Me,
                              libMesh::DenseMatrix<libMesh::Number> &Re,
                              const libMesh::Elem &element, libMesh::QGauss &q_rule,
                              const std::vector<libMesh::Real> &JxW, const std::vector<std::vector<libMesh::Real>> &phi, const std::vector<std::vector<libMesh::RealGradient>> &dphi,
                              std::vector<libMesh::dof_id_type> &dof_indices,
                              std::vector<libMesh::dof_id_type> &dof_indices_disp_x, std::vector<libMesh::dof_id_type> &dof_indices_disp_y, std::vector<libMesh::dof_id_type> &dof_indices_disp_z,
                              std::vector<libMesh::dof_id_type> &dof_indices_rot_x, std::vector<libMesh::dof_id_type> &dof_indices_rot_y, std::vector<libMesh::dof_id_type> &dof_indices_rot_z,
                              std::string element_properties_path);
void apply_boundary_condition(libMesh::EquationSystems & equation_system,
                              const std::string & system_name,
                              std::string boundary_condition_path);
void apply_load(libMesh::EquationSystems & equation_system,
                const std::string & system_name,
                std::string boundary_condition_path);
void apply_auto_spc(libMesh::EquationSystems & equation_system,
                    const std::string & system_name);
int main(int argc, char *argv[]) {
    // Command line options for vulcan
    std::string mesh_path;
    std::string boundary_condition_path;
    std::string element_properties_path;
    std::string output_path;
    for (int input = 0; input < argc; input++) {
        if (std::string(argv[input]) == "--mesh") {
            input++;
            mesh_path = argv[input];
        }
        if (std::string(argv[input]) == "--boundary_condition") {
            input++;
            boundary_condition_path = argv[input];
        }
        if (std::string(argv[input]) == "--element_properties") {
            input++;
            element_properties_path = argv[input];
        }
        if (std::string(argv[input]) == "--output") {
            input++;
            output_path = argv[input];
        }
    }
    // Initialize libMesh
    libMesh::LibMeshInit init(argc, argv);

    // Initialize mesh and read a gmsh file
    libMesh::Mesh mesh(init.comm());
    libMesh::GmshIO(mesh).read(mesh_path);
    mesh.prepare_for_use();

    // Define the elasticity system for the dynamic problem
    libMesh::EquationSystems equation_system(mesh);
    add_structural_variables(equation_system, "static");

    /*
     * Initialize the data structure for the equation system.
     */
    equation_system.init();

    assemble_system(equation_system, "static", element_properties_path);
    apply_boundary_condition(equation_system, "static", boundary_condition_path);
    apply_load(equation_system, "static", boundary_condition_path);
    apply_auto_spc(equation_system, "static");
    equation_system.get_system<libMesh::LinearImplicitSystem>("static").matrix->close();
    equation_system.get_system<libMesh::LinearImplicitSystem>("static").rhs->close();
    equation_system.get_system<libMesh::LinearImplicitSystem>("static").assemble_before_solve = false;

    equation_system.get_system<libMesh::LinearImplicitSystem>("static").solve();

    libMesh::ExodusII_IO(mesh).write_equation_systems("solution.exo", equation_system);

    std::ofstream output_file;
    output_file.open(output_path);
    for (int i = 0; i < equation_system.get_system<libMesh::LinearImplicitSystem>("static").solution->size(); i++) {
        output_file << i + 1 << ", " << equation_system.get_system<libMesh::LinearImplicitSystem>("static").solution->operator()(i) << std::endl;
    }
    output_file.close();
    return 0;
}

// Function definition
void add_structural_variables(libMesh::EquationSystems & equation_system,
                              const std::string & system_name) {
    libMesh::LinearImplicitSystem & static_system =
            equation_system.add_system<libMesh::LinearImplicitSystem>("static");
    static_system.add_vector("AUTOSPC");
    static_system.add_variable("displacement_x", libMesh::FIRST, libMesh::LAGRANGE);
    static_system.add_variable("displacement_y", libMesh::FIRST, libMesh::LAGRANGE);
    static_system.add_variable("displacement_z", libMesh::FIRST, libMesh::LAGRANGE);
    static_system.add_variable("rotation_x", libMesh::FIRST, libMesh::LAGRANGE);
    static_system.add_variable("rotation_y", libMesh::FIRST, libMesh::LAGRANGE);
    static_system.add_variable("rotation_z", libMesh::FIRST, libMesh::LAGRANGE);
}

/*
 * Constrain extra degrees of freedom required for solution.
 */
void apply_auto_spc(libMesh::EquationSystems & equation_system,
                    const std::string & libmesh_dbg_var(system_name)) {
    const libMesh::MeshBase & mesh = equation_system.get_mesh();
    libMesh::LinearImplicitSystem & static_system =
            equation_system.get_system<libMesh::LinearImplicitSystem>("static");

    for (int diag_n = 0; diag_n < static_system.get_vector("AUTOSPC").size(); diag_n++) {
        if (static_system.get_vector("AUTOSPC").operator()(diag_n) < 0.01)
            static_system.matrix->set(diag_n, diag_n, 1e15);
    }
}
/*
 * Assemble the stiffness matrix based on the information from *.elem_prop file
 */
void assemble_system(libMesh::EquationSystems & equation_system,
                     const std::string & libmesh_dbg_var(system_name),
                     std::string element_properties_path) {
    /*
     * Get the mesh connectivity, element name, and dimension
     */
    const libMesh::MeshBase & mesh = equation_system.get_mesh();
    const unsigned int dim = mesh.mesh_dimension();

    // Get a reference to the dynamic solution object
    libMesh::LinearImplicitSystem & static_system =
            equation_system.get_system<libMesh::LinearImplicitSystem>("static");

    const unsigned int disp_x_var = static_system.variable_number("displacement_x");
    const unsigned int disp_y_var = static_system.variable_number("displacement_y");
    const unsigned int disp_z_var = static_system.variable_number("displacement_z");
    const unsigned int rot_x_var = static_system.variable_number("rotation_x");
    const unsigned int rot_y_var = static_system.variable_number("rotation_y");
    const unsigned int rot_z_var = static_system.variable_number("rotation_z");
    libMesh::FEType fe_type = static_system.variable_type(disp_x_var);
    /*
     * Build a finite-element object
     */
    libMesh::UniquePtr<libMesh::FEBase> fe (libMesh::FEBase::build(dim, fe_type));

    /*
     * Define a quadrature rule for integration and attached it
     * to the finite-element object you have created.
     */
    libMesh::QGauss q_rule(dim, fe_type.default_quadrature_order());
    fe->attach_quadrature_rule(&q_rule);

    /*
     * Define references to element-specific data the I use to
     * assemble the linear system. The orientation vector is
     * what I use for discreate elements such as trusses and
     * beams. The first component if orientation vector is the
     * angle between x-z plane, and the second is the angle
     * between the projection of truss on x-z plane and x axis.
     * I use the orientation vector for rotating the stiffness
     * matrix of truss and beam in the correct position.
     */
    const std::vector<libMesh::Real> & JxW = fe->get_JxW();
    const std::vector<std::vector<libMesh::Real>> & phi = fe->get_phi();
    const std::vector<std::vector<libMesh::RealGradient>> & dphi = fe->get_dphi();

    /*
     * Get a reference to the DofMap object for this system. The
     * DofMap object handles the index translation from node and
     * element numbers to degree of freedom numbers.
     */
    const libMesh::DofMap & dof_map = static_system.get_dof_map();

    /*
     * Define data structures to contain the element stiffness, damping,
     * and mass matrix. We define another vector to store the
     * right-hand-side contribution.
     */
    libMesh::DenseMatrix<libMesh::Number> Ke;
    libMesh::DenseMatrix<libMesh::Number> Ce;
    libMesh::DenseMatrix<libMesh::Number> Me;
    libMesh::DenseMatrix<libMesh::Number> Re;
    libMesh::DenseVector<libMesh::Number> Fe;

    /*
     * The following vector, stores sum of stiffness at each degree of freedom.
     * This is used for applying auto spc command later to constrain the free
     * DOFs.
     */


    /*
     * Vector for holding degree of freedom indices for each element
     */
    std::vector<libMesh::dof_id_type> dof_indices;
    std::vector<libMesh::dof_id_type> dof_indices_disp_x;
    std::vector<libMesh::dof_id_type> dof_indices_disp_y;
    std::vector<libMesh::dof_id_type> dof_indices_disp_z;
    std::vector<libMesh::dof_id_type> dof_indices_rot_x;
    std::vector<libMesh::dof_id_type> dof_indices_rot_y;
    std::vector<libMesh::dof_id_type> dof_indices_rot_z;

    /*
     * We now loop over all the elements in the mesh and calcualte
     * element matrices and right-hand-side contributions.
     */
    libMesh::MeshBase::const_element_iterator el_itr = mesh.elements_begin();
    libMesh::MeshBase::const_element_iterator el_end = mesh.elements_end();

    /*
     * Looping over all the elements to assemble the stiffness, damping,
     * and mass matrices. I am using the element shape function and its
     * derivative along element orientation to define these matrices.
     */
    for ( ; el_itr != el_end; ++el_itr) {
        const libMesh::Elem *element = *el_itr;

        /*
         * Get the element type name. This is defined in the gmsh file. I am
         * using Nastran's element type definitions. I chose Nastran because
         * I think it is the industry standard.
         */
        std::string element_type = mesh.subdomain_name(element->subdomain_id());

        /*
         * Get the degree of freedom indices for the current element. These
         * define where in the global matrix and right-hand-side this element
         * will contribute too
         */
        dof_map.dof_indices(element, dof_indices);
        dof_map.dof_indices(element, dof_indices_disp_x, disp_x_var);
        dof_map.dof_indices(element, dof_indices_disp_y, disp_y_var);
        dof_map.dof_indices(element, dof_indices_disp_z, disp_z_var);
        dof_map.dof_indices(element, dof_indices_rot_x, rot_x_var);
        dof_map.dof_indices(element, dof_indices_rot_y, rot_y_var);
        dof_map.dof_indices(element, dof_indices_rot_z, rot_z_var);

        /*
         * Conpute the element-specific data for the current element. This
         * involves computing shape function and their derivative. We also
         * need to calculate element orientation in space for structural
         * elements such as truss and beam.
         */
        fe->reinit(element);

        /*
         * Now the calculate the FE matrices at the element level. This includes
         * calculation of stiffness, damping, and mass matrix based on the
         * element type defined in the gmsh file.
         */
        calculateElementMatrices(element_type,
                                 Ke, Ce, Me,
                                 Re,
                                 *element, q_rule,
                                 JxW, phi, dphi,
                                 dof_indices,
                                 dof_indices_disp_x, dof_indices_disp_y, dof_indices_disp_z,
                                 dof_indices_rot_x, dof_indices_rot_y, dof_indices_rot_z,
                                 element_properties_path);

        /*
         * Apply the boundary condition to the said element. This has to be done when you are assembling the system and
         * not after you are done. Changing the sparsity to assign the boundary condition can cause a mallloc error.
         * This is specified in this (https://goo.gl/959Xqa) thread. I am including one of the comments in here:
         * The error message "New nonzero at (1,9) caused a malloc!" is from PETSc.
         * PETSc is trying to do you a favor by telling you that you're allocating new
         * non-zeros in your matrix, which is slow and should be avoided.
         *
         * PETSc requires the number of nonzeros per row to be specified before the
         * matrix is assembled so that it can allocate memory efficiently. libMesh
         * handles this pre-allocation for you automatically, but libMesh assumes a
         * standard sparsity pattern based on the finite element degrees of freedom in
         * your system. In your case, you've added some extra "non-standard" non-zeros
         * into the matrix, which causes the PETSc error.
         *
         * To get around this you can either:
         *
         * 1) Not allocate the "non-standard" non-zeros in the matrix.
         *
         * 2) Tell PETSc not to report an error when a new non-zero is detected. You
         * can do this by setting MAT_NEW_NONZERO_ALLOCATION_ERR to PetscFalse by
         * calling MatSetOption on the PETSc Mat object (which you can obtain from the
         * libMesh PetscMatrix object).
         *
         * 3) Augment the libMesh sparsity pattern by calling augment_sparsity_pattern
         * on a System's dof_map. This is illustrated in libMesh example
         * miscellaneous_ex9, for example.
         */

        /*
         * Finally, I add the contribution to right places in the
         * system stiffness, damping, and mass matrices. For the
         * overal matrix, I add bunch of zeros.
         */
        static_system.matrix->add_matrix(Ke, dof_indices);

        /*
         * Check stiffness at each degree of freedom by adding the number
         * on the diagonal of Ke matrix to its corresponding DOF at SPC vector
         */
        for (short int i_diag = 0; i_diag < Ke.m(); i_diag++) {
            static_system.get_vector("AUTOSPC").add(dof_indices[i_diag], Ke.operator()(i_diag, i_diag));
        }
    }
}

/*
 * The following function calculates the element orientation and
 * the rotation matrix. The rotation matrix (Re) is used in assemble
 * stiffness function.
 */
void calculateElementOrientation(const libMesh::Elem * element,
                                 libMesh::DenseMatrix<libMesh::Real> & Re) {
    /*
     * Zero input matrix
     */
    Re.zero();
    Re.resize(2, 6);

    /*
     * Define a vector that represent the element orientation
     */
    std::vector<libMesh::Real> element_vector(3);
    element_vector[0] = element->point(1).operator()(0) - element->point(0).operator()(0);
    element_vector[1] = element->point(1).operator()(1) - element->point(0).operator()(1);
    element_vector[2] = element->point(1).operator()(2) - element->point(0).operator()(2);

    /*
     * We now calculate the direction cosines of the element. In analytic geometry,
     * the direction cosines of a vector are the cosines of the angles between the
     * vector and the three coordinate axes.
     */
    libMesh::Real element_vector_norm = sqrt(pow(element_vector[0], 2.) +
                                             pow(element_vector[1], 2.) +
                                             pow(element_vector[2], 2.));
    libMesh::Real Cx = element_vector[0] / element_vector_norm;
    libMesh::Real Cy = element_vector[1] / element_vector_norm;
    libMesh::Real Cz = element_vector[2] / element_vector_norm;

    /*
     * Now I sunstitude the direction cosine in the rotation matrix. This matrix is
     * used to rotate the 2X2 truss element stiffness matrix to 6X6 needed. Also, this
     * is used to convert the global displacements 6X1 to 2X1 needed for strain calculation
     * for the truss element. For more information on this: https://goo.gl/y5Nbj2
     * Please notice that due to the DOF numbering in libMesh, I have modified the
     * matrix defined in above source.
     */
    Re(0, 0) = Cx; Re(0, 1) = 0.; Re(0, 2) = Cy; Re(0, 3) = 0.; Re(0, 4) = Cz; Re(0, 5) = 0.;
    Re(1, 0) = 0.; Re(1, 1) = Cx; Re(1, 2) = 0.; Re(1, 3) = Cy; Re(1, 4) = 0.; Re(1, 5) = Cz;
}

/*
 * I use this function to apply the boundary condition on the
 * equation system. Please note the loads are applied through
 * a different function.
 */
void apply_boundary_condition(libMesh::EquationSystems & equation_system,
                              const std::string & system_name,
                              std::string boundary_condition_path) {
    libMesh::LinearImplicitSystem & static_system = equation_system.get_system<libMesh::LinearImplicitSystem>("static");
    libMesh::Real penalty_value = 1e15;

    /*
     * Get a constant reference to the mesh object
     */
    const libMesh::MeshBase & mesh = equation_system.get_mesh();

    /*
     * Define an iterator to go through all the nodes in the
     * boundary condition map.
     */
    Extra::boundaryCondition boundary_condition;
    boundary_condition.read(boundary_condition_path);

    std::map<unsigned int, unsigned int>::iterator bc_itr = boundary_condition.displacement_boundary.begin();
    for (bc_itr; bc_itr!=boundary_condition.displacement_boundary.end(); bc_itr++) {
        unsigned int node_number = bc_itr->first;
        /*
         * Now we need to break the DOF number to its digits. I define the
         * DOF of boundaries as an integer link 13. This means that first
         * and third degree of freedom needs to be constrained. To apply this
         * to the stiffness matrix, I need to break 13 into 1 and 3! this is
         * done by dividing 13 by 10, calculating the reminder. To calculate
         * 1, I divide 13 by 10 and cast the result into integer to remove
         * numbers after decimal. This is done until the reminder is equal
         * to zero.
         */
        unsigned int bc_dof_component;
        unsigned int bc_dof_long = bc_itr->second;
        /*
         * Continue until no digits are remaining
         */
        while (true) {
            bc_dof_component = bc_dof_long % 10;
            if (bc_dof_component > 0) {
                bc_dof_long = (int) bc_dof_long / 10;
            } else
                break;
            const libMesh::Node & current_node = mesh.node_ref(node_number);
            /*
             * Apply the dirichlet boundary condition at the corresponding
             * degrees of freedom. To get the global degree of freedom
             * associated with each node I use the .dof_number
             * function. This function takes 3 integers.
             * The first integer defines the system number (s) that I am
             * interested in. I know that my dynamic system is system 0. However,
             * I need to find a way to ask for this! The second integer (v), is the
             * number associated with the variable I am interested in. For current
             * system I have 3 variables: disp_x, disp_y, and disp_z. These are referred
             * to as '0', '1', and '2'. The third number is the component (c) of the
             * variable that I am interested in. Here, disp_x has only one component;
             * therefore, this variable only takes '0'.
             */
            static_system.matrix->set(current_node.dof_number(0, bc_dof_component - 1, 0),
                                      current_node.dof_number(0, bc_dof_component - 1, 0),
                                      penalty_value);

            /*
             * Add to SPC vector
             */
            static_system.get_vector("AUTOSPC").add(current_node.dof_number(0, bc_dof_component - 1, 0), penalty_value);
        }
    }
//    for (int i=7; i <= 11; i++) {
//        static_system.matrix->set(i, i, penalty_value);
//    }
}

/*
 * I use this function to apply load on the structure
 */
void apply_load(libMesh::EquationSystems & equation_system,
                const std::string & system_name,
                std::string boundary_condition_path) {
    libMesh::LinearImplicitSystem & static_system = equation_system.get_system<libMesh::LinearImplicitSystem>("static");

    /*
     * Get a constant referece to the mesh object
     */
    const libMesh::MeshBase & mesh = equation_system.get_mesh();

    /*
     * Define an iterator to go through all the nodes in the boundary condition map
     */
    Extra::boundaryCondition boundary_condition;

    boundary_condition.read(boundary_condition_path);
    std::map<unsigned int, std::vector<double>>::iterator bc_itr = boundary_condition.load_boundary.begin();
    for (bc_itr; bc_itr != boundary_condition.load_boundary.end(); bc_itr++) {
        unsigned int node_number = bc_itr->first;

        /*
         * Get a handle to the current node
         */
        const libMesh::Node & current_node = mesh.node_ref(node_number);
        static_system.rhs->set(current_node.dof_number(0, 0, 0), bc_itr->second[0]);
        static_system.rhs->set(current_node.dof_number(0, 1, 0), bc_itr->second[1]);
        static_system.rhs->set(current_node.dof_number(0, 2, 0), bc_itr->second[2]);
    }
}

/*
 * I use this function to calculate stiffness, damping, and mass matrix for different
 * element typre. A new element type can be added using another "if" statement.
 */
void calculateElementMatrices(std::string element_type,
                              libMesh::DenseMatrix<libMesh::Number> &Ke, libMesh::DenseMatrix<libMesh::Number> &Ce, libMesh::DenseMatrix<libMesh::Number>&Me,
                              libMesh::DenseMatrix<libMesh::Number> &Re,
                              const libMesh::Elem &element, libMesh::QGauss &q_rule,
                              const std::vector<libMesh::Real> &JxW, const std::vector<std::vector<libMesh::Real>> &phi, const std::vector<std::vector<libMesh::RealGradient>> &dphi,
                              std::vector<libMesh::dof_id_type> &dof_indices,
                              std::vector<libMesh::dof_id_type> &dof_indices_disp_x, std::vector<libMesh::dof_id_type> &dof_indices_disp_y, std::vector<libMesh::dof_id_type> &dof_indices_disp_z,
                              std::vector<libMesh::dof_id_type> &dof_indices_rot_x, std::vector<libMesh::dof_id_type> &dof_indices_rot_y, std::vector<libMesh::dof_id_type> &dof_indices_rot_z,
                              std::string element_properties_path) {

    /*
     * Read the element properteis from the *.elem_prop file
     */
    Extra::elementProperty element_property;
    element_property.read(element_properties_path);

    /*
     * Define matrices for ROD element. This is the same element as
     * two noded truss element.
     */
    if (element_type.compare("ROD") == 0) {
        /*
         * Get ROD element properties (E and A) from the input file.
         * element_property stores different properties in a map
         * data type. These can be assecced using subdomain id.
         */
        libMesh::Real elastic_modulus = element_property.material_property.at(element.subdomain_id())[0];
        libMesh::Real density = element_property.material_property.at(element.subdomain_id())[1];
        libMesh::Real cross_sectional_area = element_property.cross_section_property.at(element.subdomain_id())[0];

        /*
         * Calculate the rotation matrix and initilize system matrices.
         */
        calculateElementOrientation(&element, Re);
        Ke.resize(dof_indices_disp_x.size(), dof_indices_disp_x.size());
        Ce.resize(dof_indices_disp_x.size(), dof_indices_disp_x.size());
        Me.resize(dof_indices.size(), dof_indices.size());
        /*
         * I define different mass matrices for x, y, and z direction. The goal is
         * to calculate these and put them in the right position in the element
         * mass matrix, Me. Please look at https://goo.gl/DBXzcH for more info.
         * As mentioned in above reference, by assembeling the mass matrix in
         * this manner, there is no need to apply the rotations.
         * Me_x is the mass felt moving in the x direction, Me_y is the mass
         * felt moving in the y direction, and Me_z is the mass felt moving in
         * the z direction. As you can tell, these are all the same, but fill
         * different positions in the mass matrix.
         */
        libMesh::DenseSubMatrix<libMesh::Real> Me_x(Me), Me_y(Me), Me_z(Me);
        Me_x.reposition(0, 0, 2, 2);
        Me_y.reposition(2, 2, 2, 2);
        Me_z.reposition(4, 4, 2, 2);
        for (unsigned int qp=0; qp<q_rule.n_points(); qp++)
            for (unsigned int i=0; i < dof_indices_disp_x.size(); i++)
                for (unsigned int j=0; j < dof_indices_disp_x.size(); j++) {
                    Ke(i, j) += elastic_modulus * cross_sectional_area / element.length(0, 1) * JxW[qp] * (dphi[i][qp] * dphi[j][qp]);
                    Me_x(i, j) += density * cross_sectional_area * element.length(0, 1) * JxW[qp] * phi[i][qp] * phi[j][qp];
                    Me_y(i, j) += density * cross_sectional_area * element.length(0, 1) * JxW[qp] * phi[i][qp] * phi[j][qp];
                    Me_z(i, j) += density * cross_sectional_area * element.length(0, 1) * JxW[qp] * phi[i][qp] * phi[j][qp];
                }
        Ke.left_multiply_transpose(Re);
        Ke.right_multiply(Re);
    }
}
