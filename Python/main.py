import numpy as np

class DOMAIN():
    def __init__(self, input_file_path):
        with open(input_file_path, 'r') as input_file:
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

class DEGREE_OF_FREEDOM:
    def __init__(self):
        self.number_of_variables = 3
        self.dof = [0, 1, 2]

    def set_number_of_variables(self, number_of_variables):
        self.number_of_variables = number_of_variables
        self.dof = range(0, self.number_of_variables)

    def get_dof(self):
        self.dof = [x + self.number_of_variables * self.node_number for x in self.dof]

    def print_dof(self):
        self.get_dof()
        print self.dof

class ORIENTATION():
    def calculate_orientation(self):
        self.element_length = np.sqrt((self.node[1].x - self.node[0].x) ** 2.0 +
                                      (self.node[1].y - self.node[0].y) ** 2.0 +
                                      (self.node[1].z - self.node[0].z) ** 2.0)

        self.alpha = (self.node[1].x - self.node[0].x) / self.element_length
        self.beta = (self.node[1].y - self.node[0].y) / self.element_length
        self.gamma = (self.node[1].z - self.node[0].z) / self.element_length

class NODE(DOMAIN):
    def __init__(self):
        self.node_number = 0

    def set_coordinate(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z

class ELEMENT(DOMAIN):
    def __init__(self, node, type, material):
        

class FEM():
    def __init__(self, mesh):
        self.mesh = mesh
        for el in mesh.element:
            elmnt = dict()
            for node in mesh.element[el][1:]:
                elmnt[node] = mesh.node[node]
            element = ELEMENT(elmnt,
                              mesh.property[mesh.element[el][0]][0:2],
                              mesh.material[mesh.property[mesh.element[el][0]][2]])


mesh = DOMAIN('sample.inp')
fem = FEM(mesh)
# element = ELEMENT()

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