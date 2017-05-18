class BEAM():
    def __init__(self, element):
        self.element = element
        self.cross_sectional_area = self.element.geometry[0]
        self.modulus_of_elasticity = self.element.material[0]
        self.shear_modulus = self.element.material[1]

    def initialize(self):
        self.element.get_length()
        self.element.get_direction_cosine()
        self.get_stiffness_matrix()

    def get_stiffness_matrix(self):
        print 'beam'