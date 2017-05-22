import numpy as np

class TRUSS():
    def __init__(self, element):
        self.element = element
        self.cross_sectional_area = self.element.geometry['A']
        self.modulus_of_elasticity = self.element.material['E']

    def initialize(self):
        self.element.get_length()
        self.element.get_direction_cosine()
        self.get_stiffness_matrix()

    def get_stiffness_matrix(self):
        self.element.matrix['stiffness'] = self.modulus_of_elasticity * self.cross_sectional_area / self.element.length * np.matrix([[1, -1], [-1, 1]])
        rotation_matrix = np.matrix([[self.element.alpha, self.element.beta, self.element.gamma, 0., 0., 0., 0., 0., 0., 0., 0., 0.],
                                     [0., 0., 0., 0., 0., 0., self.element.alpha, self.element.beta, self.element.gamma, 0., 0., 0.]])
        self.element.matrix['stiffness'] = rotation_matrix.T * self.element.matrix['stiffness'] * rotation_matrix