import numpy as np

class DEGREE_OF_FREEDOM:
    def __init__(self):
        self.dof = [0, 1, 2]
        self.number_of_variables = 3

    def set_number_of_variables(self, number_of_variables):
        self.number_of_variables = number_of_variables
        self.dof = range(0, self.number_of_variables)

    def get_dof(self):
        self.dof = [x + self.number_of_variables * self.node_number for x in self.dof]

    def print_dof(self):
        self.get_dof()
        print self.dof

class NODE(DEGREE_OF_FREEDOM):
    def __init__(self):
        DEGREE_OF_FREEDOM.__init__(self)
        self.node_number = 0

    def set_coordinate(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z

class ELEMENT():
    def __init__(self, node_list, quadrature_rule=1):
        self.node = node_list
        self.number_of_nodes = len(node_list)
        self.quadrature_rule = quadrature_rule
        self.calculate_jacobian()
        self.get_quadrature_points()
        self.get_shape_functions()
        self.get_shape_functions_derivative()
        self.get_shape_function_at_quadrature_points()
        self.get_shape_function_derivative_at_quadrature_points()
        self.calculate_stiffness_matrix()

    def calculate_jacobian(self):
        self.jacobian = np.sqrt(((self.node[1].x - self.node[0].x) / 2.)**2. + \
                                ((self.node[1].y - self.node[0].y) / 2.)**2. + \
                                ((self.node[1].z - self.node[0].z) / 2.)**2.)

    def get_quadrature_points(self):
        if (self.number_of_nodes == 2):
            self.quadrature_points = [-np.sqrt(1. / 3.), np.sqrt(1. / 3.)]
            self.quadrature_points_weight = [1., 1.]

    def get_shape_functions(self):
        if (self.number_of_nodes == 2):
            self.shape_functions = [lambda zeta: (1. - zeta) / 2.,
                                    lambda zeta: (1. + zeta) / 2.]

    def get_shape_functions_derivative(self):
        if (self.number_of_nodes == 2):
            self.shape_functions_derivative = [lambda zeta: -1. / 2.,
                                               lambda zeta: 1. / 2.]

    def get_shape_function_at_quadrature_points(self):
        if (self.number_of_nodes == 2):
            self.shape_functions_at_quadrature_points = np.matrix([[self.shape_functions[0](self.quadrature_points[0]),
                                                                    self.shape_functions[1](self.quadrature_points[0])],
                                                                   [self.shape_functions[0](self.quadrature_points[1]),
                                                                    self.shape_functions[1](self.quadrature_points[1])]])

    def get_shape_function_derivative_at_quadrature_points(self):
        if (self.number_of_nodes == 2):
            self.shape_functions_derivative_at_quadrature_points = np.matrix([[self.shape_functions_derivative[0](self.quadrature_points[0]),
                                                                               self.shape_functions_derivative[1](self.quadrature_points[0])],
                                                                              [self.shape_functions_derivative[0](self.quadrature_points[1]),
                                                                               self.shape_functions_derivative[1](self.quadrature_points[1])]])

    def calculate_stiffness_matrix(self):
        self.stiffness_matrix = np.zeros([2, 2])
        for i_quadrature_points in range(0, len(self.quadrature_points)):
            for i_shape_function in range(0, len(self.shape_functions)):
                for j_shape_function in range(0, len(self.shape_functions)):
                    self.stiffness_matrix[i_shape_function, j_shape_function] += \
                                                                                self.shape_functions_derivative_at_quadrature_points[i_quadrature_points, i_shape_function] * \
                                                                                self.shape_functions_derivative_at_quadrature_points[i_quadrature_points, j_shape_function] * \
                                                                                1. / self.jacobian * \
                                                                                self.quadrature_points_weight[i_quadrature_points]



node = []
node.append(NODE())
node.append(NODE())

node[0].set_coordinate(0, 0, 0)
node[0].node_number = 0

node[1].set_coordinate(1, 0, 0)
node[1].node_number = 1

element = ELEMENT(node)
print element.stiffness_matrix

