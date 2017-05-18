from NODE import NODE
from TRUSS import TRUSS
from BEAM import BEAM

import numpy as np

class ELEMENT(NODE):
    def __init__(self):
        self.id = None
        self.node = []
        self.type = None
        self.geometry = None
        self.material = None
        self.stiffness_matrix = None

    def attach_nodes(self, nodes):
        try:
            for node in nodes:
                self.node.append(node)
        except TypeError:
            self.node.append(nodes)

    def get_length(self):
        self.length = np.sqrt((self.node[-1].x - self.node[0].x) ** 2. +
                              (self.node[-1].y - self.node[0].y) ** 2. +
                              (self.node[-1].z - self.node[0].z) ** 2.)

    def get_direction_cosine(self):
        self.get_length()
        self.alpha = (self.node[-1].x - self.node[0].x) / self.length
        self.beta = (self.node[-1].y - self.node[0].y) / self.length
        self.gamma = (self.node[-1].z - self.node[0].z) / self.length

    def initialize(self):
        if self.type == 'TRUSS':
            truss = TRUSS(self)
            truss.initialize()
        elif self.type == 'BEAM':
            beam = BEAM(self)
            beam.initialize()

    # def get_stiffness_matrix(self):
    #     if self.type == 'TRUSS':
    #         truss = TRUSS(self)
    #         truss.get_stiffness_matrix()
    #     elif self.type == 'BEAM':
    #         beam = BEAM(self)
    #         beam.get_stiffness_matrix()