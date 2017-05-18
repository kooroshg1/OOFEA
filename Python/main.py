# Add src to path
import sys
sys.path.insert(0, 'src/mesh')
sys.path.insert(0, 'src/base')
import MESH


import numpy as np
import scipy.sparse as sps
import scipy.sparse.linalg as spsol

class NODE_BASE():
    def get_nodal_degree_of_freedom(self):
        self.nodal_degree_of_freedom = []
        for node_number in self.node:
            for dof in range(0, self.dimension):
                self.nodal_degree_of_freedom.append(self.dimension * node_number + dof)

class ELEMENT_1D_BASE(NODE_BASE):
    def calculate_length(self):
        self.length = np.sqrt((self.node[1][0] - self.node[0][0]) ** 2. +
                              (self.node[1][1] - self.node[0][1]) ** 2. +
                              (self.node[1][2] - self.node[0][2]) ** 2.)

    def calculate_direction_cosines(self):
        self.alpha = (self.node[1][0] - self.node[0][0]) / self.length
        self.beta = (self.node[1][1] - self.node[0][1]) / self.length
        self.gamma = (self.node[1][2] - self.node[0][2]) / self.length

    def calculate_transformation_matrix(self):
        node[1][0] = self.node[1][0] - self.node[0][0]
        node[1][1] = self.node[1][1] - self.node[0][1]
        node[1][2] = self.node[1][2] - self.node[0][2]
        node[0] = [0., 0., 0.]

        L = np.sqrt((node[1][0] - node[0][0]) ** 2. +
                    (node[1][1] - node[0][1]) ** 2. +
                    (node[1][2] - node[0][2]) ** 2.)
        x12 = np.sqrt(node[1][0] ** 2. + node[1][1] ** 2.)
        x11 = np.sqrt(node[0][0] ** 2. + node[0][1] ** 2.)

        cosAlpha = (node[1][0] - node[0][0]) / (x12 - x11)
        sinAlpha = (node[1][1] - node[0][1]) / (x12 - x11)
        cosBeta = (x12 - x11) / L
        sinBeta = (node[1][2] - node[0][2]) / L

        LAMBDA = np.matrix([[cosAlpha * cosBeta, sinAlpha, cosAlpha * sinBeta],
                            [-sinAlpha * cosBeta, cosAlpha, -sinAlpha * sinBeta],
                            [-sinBeta, 0., cosBeta]])

        self.transformation_matrix = np.kron(np.eye(4, 4), LAMBDA)

class TRUSS(ELEMENT_1D_BASE):
    def __init__(self, node, property_list):
        self.node = node
        self.type = element_type
        self.property_list = property_list
        self.calculate_length()
        self.calculate_direction_cosines()
        self.calculate_stiffness()

    def calculate_stiffness(self):
        self.stiffness_matrix = np.matrix([[1, -1], [-1, 1]])
        self.stiffness_matrix = self.property_list['YOUNG'] * self.property_list[
            'AREA'] / self.length * self.stiffness_matrix
        R = np.matrix([[self.alpha, self.beta, self.gamma, 0., 0., 0., 0., 0., 0., 0., 0., 0.],
                       [0., 0., 0., 0., 0., 0., self.alpha, self.beta, self.gamma, 0., 0., 0.]])
        self.stiffness_matrix = R.T * self.stiffness_matrix * R

class BEAM(ELEMENT_1D_BASE):
    def __init__(self, node, properties):
        self.node = node
        self.properties = properties
        self.calculate_length()
        self.calculate_transformation_matrix()
        self.calculate_local_element_stiffness_matrix()
        self.calculate_global_element_stiffness_matrix()

    def calculate_local_element_stiffness_matrix(self):
        self.stiffness_matrix = np.zeros([12, 12])
        self.E = self.properties['YOUNG']
        self.G = self.properties['G']
        self.area = self.properties['AREA']
        self.Iy = self.properties['Iy']
        self.Iz = self.properties['Iz']
        self.J = self.properties['J']
        self.a = self.length / 2.

        self.stiffness_matrix[0, 0] = self.area * self.E / (2. * self.a)
        self.stiffness_matrix[0, 6] = -self.area * self.E / (2. * self.a)

        self.stiffness_matrix[1, 1] = 3 * self.E * self.Iz / (2. * self.a ** 3.)
        self.stiffness_matrix[1, 5] = 3 * self.E * self.Iz / (2. * self.a ** 2.)
        self.stiffness_matrix[1, 7] = -3 * self.E * self.Iz / (2. * self.a ** 3.)
        self.stiffness_matrix[1, 11] = 3 * self.E * self.Iz / (2. * self.a ** 2.)

        self.stiffness_matrix[2, 2] = 3 * self.E * self.Iy / (2. * self.a ** 3.)
        self.stiffness_matrix[2, 4] = -3 * self.E * self.Iy / (2. * self.a ** 2.)
        self.stiffness_matrix[2, 8] = -3 * self.E * self.Iy / (2. * self.a ** 3.)
        self.stiffness_matrix[2, 10] = -3 * self.E * self.Iy / (2. * self.a ** 2.)

        self.stiffness_matrix[3, 3] = self.G * self.J / (2. * self.a)
        self.stiffness_matrix[3, 9] = -self.G * self.J / (2. * self.a)

        self.stiffness_matrix[4, 4] = 2. * self.E * self.Iy / self.a
        self.stiffness_matrix[4, 8] = 3. * self.E * self.Iy / (2. * self.a ** 2.)
        self.stiffness_matrix[4, 10] = self.E * self.Iy / self.a

        self.stiffness_matrix[5, 5] = 2. * self.E * self.Iz / self.a
        self.stiffness_matrix[5, 7] = -3. * self.E * self.Iz / (2. * self.a ** 2.)
        self.stiffness_matrix[5, 11] = self.E * self.Iz / self.a

        self.stiffness_matrix[6, 6] = self.area * self.E / (2. * self.a)

        self.stiffness_matrix[7, 7] = 3 * self.E * self.Iz / (2. * self.a ** 3.)
        self.stiffness_matrix[7, 11] = -3 * self.E * self.Iz / (2. * self.a ** 2.)

        self.stiffness_matrix[8, 8] = 3 * self.E * self.Iy / (2. * self.a ** 3.)
        self.stiffness_matrix[8, 10] = 3 * self.E * self.Iy / (2. * self.a ** 2.)

        self.stiffness_matrix[9, 9] = self.G * self.J / (2. * self.a ** 3.)

        self.stiffness_matrix[10, 10] = 2 * self.E * self.Iy / self.a

        self.stiffness_matrix[11, 11] = 2 * self.E * self.Iz / self.a

        self.stiffness_matrix = self.stiffness_matrix + self.stiffness_matrix.T - np.diag(self.stiffness_matrix.diagonal())

    def calculate_global_element_stiffness_matrix(self):
        self.stiffness_matrix = self.transformation_matrix * self.stiffness_matrix

class ELEMENT(TRUSS, BEAM):
    def __init__(self, node, element_type, property_list):
        self.dimension = 6
        self.node = node
        self.type = element_type
        self.property_list = property_list
        self.calculate_stiffness()
        self.get_nodal_degree_of_freedom()

    def calculate_stiffness(self):
        if self.type == 'TRUSS':
            truss = TRUSS(self.node, self.property_list)
            self.stiffness_matrix = truss.stiffness_matrix
        elif self.type == 'BEAM':
            beam = BEAM(self.node, self.property_list)
            self.stiffness_matrix = beam.stiffness_matrix

class BOUNDARY_CONDITION_BASE():
    def __init__(self):
        self.node_set = None
        self.disp_x = None
        self.disp_y = None
        self.disp_z = None
        self.rot_x = None
        self.rot_y = None
        self.rot_z = None
        self.F_x = None
        self.F_y = None
        self.F_z = None
        self.M_x = None
        self.M_y = None
        self.M_z = None
        self.AUTOSPC = None

class BOUNDARY_CONDITION(BOUNDARY_CONDITION_BASE):
    def add_displacement_boundary_condition(self, node_set, disp_x, disp_y, disp_z):
        self.node_set = node_set
        self.disp_x = disp_x
        self.disp_y = disp_y
        self.disp_z = disp_z

    def add_rotation_boundary_condition(self, node_set, rot_x, rot_y, rot_z):
        self.node_set = node_set
        self.rot_x = rot_x
        self.rot_y = rot_y
        self.rot_z = rot_z

    def add_force_boundary_condition(self, node_set, F_x, F_y, F_z):
        self.node_set = node_set
        self.F_x = F_x
        self.F_y = F_y
        self.F_z = F_z

    def add_moment_boundary_condition(self, node_set, M_x, M_y, M_z):
        self.node_set = node_set
        self.M_x = M_x
        self.M_y = M_y
        self.M_z = M_z

    def apply_displacement_boundary_condition(self):
        for bc in mesh.displacement_boundary_condition:
            for node in mesh.node_set[bc.node_set]:
                dof = [0, 1, 2]
                self.dof = [i + node * 6 for i in dof]
                self.colomns = [0 for i in dof]
                self.penalty = [1e20, 1e20, 1e20]
                self.displacement_value = [bc.disp_x, bc.disp_y, bc.disp_z]
                self.stiffness_matrix += sps.coo_matrix((self.penalty, (self.dof, self.dof)), shape=self.stiffness_matrix.shape)
                self.rhs += sps.coo_matrix((self.displacement_value, (self.dof, self.colomns)), shape=self.rhs.shape)

    def apply_rotation_boundary_condition(self):
        for bc in mesh.rotation_boundary_condition:
            for node in mesh.node_set[bc.node_set]:
                dof = [3, 4, 5]
                self.dof = [i + node * 6 for i in dof]
                self.colomns = [0 for i in dof]
                self.penalty = [1e20, 1e20, 1e20]
                self.rotation_value = [bc.rot_x, bc.rot_y, bc.rot_z]
                self.stiffness_matrix += sps.coo_matrix((self.penalty, (self.dof, self.dof)), shape=self.stiffness_matrix.shape)
                self.rhs += sps.coo_matrix((self.rotation_value, (self.dof, self.colomns)), shape=self.rhs.shape)

    def apply_force_boundary_condition(self):
        for bc in mesh.force_boundary_condition:
            for node in mesh.node_set[bc.node_set]:
                dof = [0, 1, 2]
                self.dof = [i + node * 6 for i in dof]
                self.colomns = [0 for i in dof]
                self.force_value = [bc.F_x, bc.F_y, bc.F_z]
                self.rhs += sps.coo_matrix((self.force_value, (self.dof, self.colomns)), shape=self.rhs.shape)

    def apply_moment_boundary_condition(self):
        for bc in mesh.moment_boundary_condition:
            for node in mesh.node_set[bc.node_set]:
                dof = [0, 1, 2]
                self.dof = [i + node * 6 for i in dof]
                self.colomns = [0 for i in dof]
                self.moment_value = [bc.M_x, bc.M_y, bc.M_z]
                self.rhs += sps.coo_matrix((self.moment_value, (self.dof, self.colomns)), shape=self.rhs.shape)

    def apply_auto_spc(self):
        if self.mesh.AUTOSPC:
            for diag in range(0, self.stiffness_matrix.shape[0]):
                if self.stiffness_matrix[diag, diag] < 0.1:
                    self.stiffness_matrix[diag, diag] = 1e20

class MODEL(BOUNDARY_CONDITION):
    def __init__(self, mesh):
        self.mesh = mesh
        self.stiffness_matrix = sps.coo_matrix((len(self.mesh.node) * self.mesh.dimension, len(self.mesh.node) * self.mesh.dimension))
        self.rhs = sps.coo_matrix((len(self.mesh.node) * self.mesh.dimension, 1))

    def add_matrix(self, matrix, dof):
        combinations = [(x, y) for x in dof for y in dof]
        rows = [x for (x, y) in combinations]
        columns = [y for (x, y) in combinations]
        matrix = matrix.reshape(-1, 1)
        matrix = [float(x) for x in matrix]
        self.stiffness_matrix += sps.coo_matrix((matrix, (rows, columns)), shape=self.stiffness_matrix.shape)

    def add_rhs(self, vector, dof):
        columns = [0 for y in dof]
        rows = [x for x in dof]
        vector = [float(x) for x in vector]
        self.rhs += sps.coo_matrix((vector, (rows, columns)), shape=self.rhs.shape)

    def solve(self):
        self.solution = spsol.spsolve(self.stiffness_matrix, self.rhs)

mesh = MESH.MESH()
mesh.mesh_file_path = 'sample.inp'
mesh.read_mesh()

# fea = MODEL(mesh)
# # bc = BOUNDARY_CONDITION(mesh)
# # bc.apply_displacement_boundary_condition()
# # fea.apply_displacement_boundary_condition()
#
# output = open('matrix.txt', 'w')
# for el in mesh.element:
#     property_id = mesh.element[el][0]
#     material_id = mesh.property[property_id][1]
#     element_type = mesh.property[mesh.element[el][0]][0]
#
#     node = dict()
#     property_list = dict()
#
#     for nd in mesh.element[el][1:]:
#         node[nd] = mesh.node[nd]
#
#     if element_type == 'TRUSS':
#         property_list['YOUNG'] = mesh.material[material_id][0]
#         property_list['DENSITY'] = mesh.material[material_id][1]
#         property_list['AREA'] = mesh.property[property_id][2]
#     elif element_type == 'BEAM':
#         property_list['YOUNG'] = mesh.material[material_id][0]
#         property_list['G'] = mesh.material[material_id][1]
#         property_list['DENSITY'] = mesh.material[material_id][2]
#         property_list['AREA'] = mesh.property[property_id][2]
#         property_list['J'] = mesh.property[property_id][3]
#         property_list['Ix'] = mesh.property[property_id][4]
#         property_list['Iy'] = mesh.property[property_id][5]
#         property_list['Iz'] = mesh.property[property_id][6]
#
#     element = ELEMENT(node, element_type, property_list)
#
#     fea.add_matrix(element.stiffness_matrix, element.nodal_degree_of_freedom)
#
#     np.savetxt(output, element.stiffness_matrix, fmt='%-12.4E')
#     output.write('\n')
#
# fea.apply_displacement_boundary_condition()
# fea.apply_rotation_boundary_condition()
# fea.apply_force_boundary_condition()
# fea.apply_auto_spc()
# np.savetxt(output, fea.stiffness_matrix.todense(), fmt='%-12.4E')
# output.close()
#
# # print fea.stiffness_matrix
# # print fea.rhs.todense()
#
# fea.solve()
# print fea.solution
#
#
# class vehicle():
#     def __init__(self):
#         self.value = None
#         self.is_second = None
#         self.sold = None
#
#     def mark_as_sold(self):
#         self.sold = True
#
#     def is_sold(self):
#         print self.sold
#
#     def get_value(self):
#         print self.value
#
# class car(vehicle):
#     def __init__(self, value, is_second, sold):
#         self.value = value
#         self.is_second = is_second
#         self.sold = sold
