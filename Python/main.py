import numpy as np
import scipy.sparse as sps

class DOMAIN():
    def __init__(self, input_file_path):
        self.input_file_path = input_file_path
        self.read_mesh()

    def read_mesh(self):
        self.dimension = 3
        with open(self.input_file_path, 'r') as input_file:
            self.node = dict()
            self.element = dict()
            self.property = dict()
            self.material = dict()
            for line in input_file:
                if line[0:4] == 'NODE':
                    line_split = line.split(',')
                    self.node[int(line_split[1])] = []
                    for line_split_element in line_split[2:]:
                        self.node[int(line_split[1])].append(float(line_split_element))
                elif line[0:7] == 'ELEMENT':
                    line_split = line.split(',')
                    self.element[int(line_split[1])] = []
                    for line_split_element in line_split[2:]:
                        self.element[int(line_split[1])].append(int(line_split_element))
                elif line[0:8] == 'PROPERTY':
                    line_split = line.split(',')
                    self.property[int(line_split[1])] = []
                    for line_split_element in line_split[2:]:
                        try:
                            self.property[int(line_split[1])].append(float(line_split_element))
                        except ValueError:
                            self.property[int(line_split[1])].append(line_split_element)
                elif line[0:8] == 'MATERIAL':
                    line_split = line.split(',')
                    self.material[int(line_split[1])] = []
                    for line_split_element in line_split[2:]:
                        self.material[int(line_split[1])].append(float(line_split_element))

class ELEMENT_BASE():
    def __init__(self):
        self.dim = None
        self.node_list = None
        self.type = None
        self.property_list = None
        self.nodal_degree_of_freedom = None

    def get_nodal_degree_of_freedom(self):
        self.nodal_degree_of_freedom = []
        for node_number in self.node_list:
            for dof in range(0, self.dim):
                self.nodal_degree_of_freedom.append(self.dim * node_number + dof)

class ELEMENT_1D(ELEMENT_BASE):
    def __init__(self):
        self.dim = None
        self.node_list = None
        self.type = None
        self.property_list = None

    def calculate_length(self):
        self.length = np.sqrt((self.node_list[self.node_list.keys()[0]][0] - self.node_list[self.node_list.keys()[1]][0])**2.0 +
                              (self.node_list[self.node_list.keys()[0]][1] - self.node_list[self.node_list.keys()[1]][1])**2.0 +
                              (self.node_list[self.node_list.keys()[0]][2] - self.node_list[self.node_list.keys()[1]][2])**2.0)

    def calculate_orientation(self):
        self.alpha = (self.node_list[self.node_list.keys()[1]][0] - self.node_list[self.node_list.keys()[0]][0]) / self.length
        self.beta = (self.node_list[self.node_list.keys()[1]][1] - self.node_list[self.node_list.keys()[0]][1]) / self.length
        self.gamma = (self.node_list[self.node_list.keys()[1]][2] - self.node_list[self.node_list.keys()[0]][2]) / self.length

class ELEMENT(ELEMENT_1D):
    def __init__(self, node_list, element_type, property_list):
        self.dim = 3
        self.node_list = node_list
        self.type = element_type
        self.property_list = property_list
        self.calculate_length()
        self.calculate_orientation()

    def calculate_stiffness(self):
        if self.type == 'TRUSS':
            self.stiffness_matrix = np.matrix([[1, -1],[-1, 1]])
            self.stiffness_matrix = self.property_list['YOUNG'] * self.property_list['CROSS_A'] / self.length * self.stiffness_matrix
            R = np.matrix([[self.alpha, self.beta, self.gamma, 0., 0., 0.],
                           [0., 0., 0., self.alpha, self.beta, self.gamma]])
            self.stiffness_matrix = R.T * self.stiffness_matrix * R

class MODEL():
    def __init__(self, mesh):
        self.mesh = mesh
        self.stiffness_matrix = sps.coo_matrix((len(self.mesh.node) * self.mesh.dimension, len(self.mesh.node) * self.mesh.dimension))

    def add_matrix(self, matrix, dof):
        combinations = [(x, y) for x in dof for y in dof]
        rows = [x for (x, y) in combinations]
        columns = [y for (x, y) in combinations]
        matrix = matrix.reshape(-1, 1)
        matrix = [float(x) for x in matrix]
        self.stiffness_matrix += sps.coo_matrix((matrix, (rows, columns)), shape=self.stiffness_matrix.shape)


mesh = DOMAIN('sample.inp')
fea = MODEL(mesh)

# data = np.matrix([[1, 2], [3, 1]])
# rows = [1, 2, 3, 1]
# columns = [1, 2, 3, 4]
# data = data.reshape(-1, 1)
# data =
# print data
# mat = sps.coo_matrix((data, (rows, columns)))

output = open('matrix.txt', 'w')
for el in mesh.element:
    property_id = mesh.element[el][0]
    material_id = mesh.property[property_id][2]

    element_type = mesh.property[mesh.element[el][0]][0]

    node_list = dict()
    property_list = dict()

    for node in mesh.element[el][1:]:
        node_list[node] = mesh.node[node]

    if element_type == 'TRUSS':
        property_list['YOUNG'] = mesh.material[material_id][0]
        property_list['DENSITY'] = mesh.property[material_id][1]
        property_list['CROSS_A'] = mesh.property[property_id][1]

    element = ELEMENT(node_list, element_type, property_list)
    element.calculate_stiffness()
    element.get_nodal_degree_of_freedom()

    fea.add_matrix(element.stiffness_matrix, element.nodal_degree_of_freedom)

    np.savetxt(output, element.stiffness_matrix, fmt='%-12.4E')
    output.write('\n')

np.savetxt(output, fea.stiffness_matrix.todense(), fmt='%12.4E')
output.close()


class vehicle():
    def __init__(self):
        self.value = None
        self.is_second = None
        self.sold = None

    def mark_as_sold(self):
        self.sold = True

    def is_sold(self):
        print self.sold

    def get_value(self):
        print self.value

class car(vehicle):
    def __init__(self, value, is_second, sold):
        self.value = value
        self.is_second = is_second
        self.sold = sold

# class ELEMENT(ORIENTATION):
#     def __init__(self, node_list, quadrature_rule=1):
#         self.node = node_list
#         self.number_of_nodes = len(node_list)
#         self.quadrature_rule = quadrature_rule
#         self.calculate_jacobian()
#         self.get_quadrature_points()
#         self.get_shape_functions()
#         self.get_shape_functions_derivative()
#         self.get_shape_function_at_quadrature_points()
#         self.get_shape_function_derivative_at_quadrature_points()
#         self.calculate_local_element_stiffness_matrix()
#         self.calculate_orientation()
#         self.calculate_global_element_stiffness_matrix()
#
#
#     def calculate_jacobian(self):
#         self.jacobian = np.sqrt(((self.node[1].x - self.node[0].x) / 2.)**2. + \
#                                 ((self.node[1].y - self.node[0].y) / 2.)**2. + \
#                                 ((self.node[1].z - self.node[0].z) / 2.)**2.)
#
#     def get_quadrature_points(self):
#         if (self.number_of_nodes == 2):
#             self.quadrature_points = [-np.sqrt(1. / 3.), np.sqrt(1. / 3.)]
#             self.quadrature_points_weight = [1., 1.]
#
#     def get_shape_functions(self):
#         if (self.number_of_nodes == 2):
#             self.shape_functions = [lambda zeta: (1. - zeta) / 2.,
#                                     lambda zeta: (1. + zeta) / 2.]
#
#     def get_shape_functions_derivative(self):
#         if (self.number_of_nodes == 2):
#             self.shape_functions_derivative = [lambda zeta: -1. / 2.,
#                                                lambda zeta: 1. / 2.]
#
#     def get_shape_function_at_quadrature_points(self):
#         if (self.number_of_nodes == 2):
#             self.shape_functions_at_quadrature_points = np.matrix([[self.shape_functions[0](self.quadrature_points[0]),
#                                                                     self.shape_functions[1](self.quadrature_points[0])],
#                                                                    [self.shape_functions[0](self.quadrature_points[1]),
#                                                                     self.shape_functions[1](self.quadrature_points[1])]])
#
#     def get_shape_function_derivative_at_quadrature_points(self):
#         if (self.number_of_nodes == 2):
#             self.shape_functions_derivative_at_quadrature_points = np.matrix([[self.shape_functions_derivative[0](self.quadrature_points[0]),
#                                                                                self.shape_functions_derivative[1](self.quadrature_points[0])],
#                                                                               [self.shape_functions_derivative[0](self.quadrature_points[1]),
#                                                                                self.shape_functions_derivative[1](self.quadrature_points[1])]])
#
#     def calculate_local_element_stiffness_matrix(self):
#         self.local_element_stiffness_matrix = np.zeros([2, 2])
#         for i_quadrature_points in range(0, len(self.quadrature_points)):
#             for i_shape_function in range(0, len(self.shape_functions)):
#                 for j_shape_function in range(0, len(self.shape_functions)):
#                     self.local_element_stiffness_matrix[i_shape_function, j_shape_function] += \
#                                                                                                self.shape_functions_derivative_at_quadrature_points[i_quadrature_points, i_shape_function] * \
#                                                                                                self.shape_functions_derivative_at_quadrature_points[i_quadrature_points, j_shape_function] * \
#                                                                                                1. / self.jacobian * \
#                                                                                                self.quadrature_points_weight[i_quadrature_points]
#
#     def calculate_global_element_stiffness_matrix(self):
#         rotation_matrix = np.matrix([[self.alpha, self.beta, self.gamma, 0., 0., 0.],
#                                      [0., 0., 0., self.alpha, self.beta, self.gamma]])
#         self.global_element_stiffness_matrix = rotation_matrix.T * \
#                                                self.local_element_stiffness_matrix * \
#                                                rotation_matrix