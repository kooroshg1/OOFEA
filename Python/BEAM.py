import numpy as np

class BEAM():
    def __init__(self, element):
        self.element = element
        self.cross_sectional_area = self.element.geometry['A']
        self.modulus_of_elasticity = self.element.material['E']
        self.shear_modulus = self.element.material['G']
        self.Ix = self.element.geometry['Ix']
        self.Iy = self.element.geometry['Iy']
        self.Iz = self.element.geometry['Iz']
        self.J = self.element.geometry['J']

    def initialize(self):
        self.element.get_length()
        self.element.get_direction_cosine()
        self.element.calculate_transformation_matrix()
        self.get_stiffness_matrix()

    def get_stiffness_matrix(self):
        self.stiffness_matrix = np.zeros([12, 12])
        self.a = self.element.length / 2.

        self.stiffness_matrix[0, 0] = self.cross_sectional_area * self.modulus_of_elasticity / (2. * self.a)
        self.stiffness_matrix[0, 6] = -self.cross_sectional_area * self.modulus_of_elasticity / (2. * self.a)

        self.stiffness_matrix[1, 1] = 3 * self.modulus_of_elasticity * self.Iz / (2. * self.a ** 3.)
        self.stiffness_matrix[1, 5] = 3 * self.modulus_of_elasticity * self.Iz / (2. * self.a ** 2.)
        self.stiffness_matrix[1, 7] = -3 * self.modulus_of_elasticity * self.Iz / (2. * self.a ** 3.)
        self.stiffness_matrix[1, 11] = 3 * self.modulus_of_elasticity * self.Iz / (2. * self.a ** 2.)

        self.stiffness_matrix[2, 2] = 3 * self.modulus_of_elasticity * self.Iy / (2. * self.a ** 3.)
        self.stiffness_matrix[2, 4] = -3 * self.modulus_of_elasticity * self.Iy / (2. * self.a ** 2.)
        self.stiffness_matrix[2, 8] = -3 * self.modulus_of_elasticity * self.Iy / (2. * self.a ** 3.)
        self.stiffness_matrix[2, 10] = -3 * self.modulus_of_elasticity * self.Iy / (2. * self.a ** 2.)

        self.stiffness_matrix[3, 3] = self.shear_modulus * self.J / (2. * self.a)
        self.stiffness_matrix[3, 9] = -self.shear_modulus * self.J / (2. * self.a)

        self.stiffness_matrix[4, 4] = 2. * self.modulus_of_elasticity * self.Iy / self.a
        self.stiffness_matrix[4, 8] = 3. * self.modulus_of_elasticity * self.Iy / (2. * self.a ** 2.)
        self.stiffness_matrix[4, 10] = self.modulus_of_elasticity * self.Iy / self.a

        self.stiffness_matrix[5, 5] = 2. * self.modulus_of_elasticity * self.Iz / self.a
        self.stiffness_matrix[5, 7] = -3. * self.modulus_of_elasticity * self.Iz / (2. * self.a ** 2.)
        self.stiffness_matrix[5, 11] = self.modulus_of_elasticity * self.Iz / self.a

        self.stiffness_matrix[6, 6] = self.cross_sectional_area * self.modulus_of_elasticity / (2. * self.a)

        self.stiffness_matrix[7, 7] = 3 * self.modulus_of_elasticity * self.Iz / (2. * self.a ** 3.)
        self.stiffness_matrix[7, 11] = -3 * self.modulus_of_elasticity * self.Iz / (2. * self.a ** 2.)

        self.stiffness_matrix[8, 8] = 3 * self.modulus_of_elasticity * self.Iy / (2. * self.a ** 3.)
        self.stiffness_matrix[8, 10] = 3 * self.modulus_of_elasticity * self.Iy / (2. * self.a ** 2.)

        self.stiffness_matrix[9, 9] = self.shear_modulus * self.J / (2. * self.a ** 3.)

        self.stiffness_matrix[10, 10] = 2 * self.modulus_of_elasticity * self.Iy / self.a

        self.stiffness_matrix[11, 11] = 2 * self.modulus_of_elasticity * self.Iz / self.a

        self.stiffness_matrix = self.stiffness_matrix + self.stiffness_matrix.T - np.diag(self.stiffness_matrix.diagonal())

        self.element.stiffness_matrix = self.element.transformation_matrix * self.stiffness_matrix