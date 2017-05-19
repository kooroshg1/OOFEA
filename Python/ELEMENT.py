from NODE import NODE
from TRUSS import TRUSS
from BEAM import BEAM

import numpy as np

class ELEMENT(NODE):
    def __init__(self):
        self.id = None
        self.node = []
        self.type = None
        self.geometry = dict()
        self.material = dict()
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

    def calculate_transformation_matrix(self):
        # Source: https://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d
        self.get_length()
        orientation = np.array([self.node[-1].x - self.node[0].x,
                                self.node[-1].y - self.node[0].y,
                                self.node[-1].z - self.node[0].z]) / self.length
        xorientation = np.array([1, 0, 0])
        v = np.cross(orientation, xorientation)
        s = np.linalg.norm(v)
        c = np.dot(orientation, xorientation)
        vx = np.matrix([[0., -v[2], v[1]],
                        [v[2], 0., -v[0]],
                        [-v[1], v[0], 0.]])
        if c == -1:
            self.transformation_matrix = np.matrix([[-1, 0, 0], [0, -1, 0], [0, 0, 1]])
        else:
            self.transformation_matrix = np.eye(3) + vx + vx ** 2 * 1. / (1. + c)

        self.transformation_matrix = np.kron(np.eye(4, 4), self.transformation_matrix)

