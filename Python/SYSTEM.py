import scipy.sparse as sps
from sets import Set

class LinearElasticity():
    def __init__(self):
        self.matrix = dict()
        self.vector = dict()
        self.elements = None
        return None

    def add_elements(self, elements):
        self.elements = elements

    def add_matrix(self, matrixName):
        self.matrix[matrixName] = sps.coo_matrix((24, 24))

    def add_vector(self, vectorName):
        self.vector[vectorName] = sps.coo_matrix((24, 1))

    def set_matrix(self, systemMatrixName, elementMatrixName):
        for element in self.elements:
            dof = []
            for node in element.node:
                dof.extend(node.dof)
            combinations = [(x, y) for x in dof for y in dof]
            rows = [x for (x, y) in combinations]
            columns = [y for (x, y) in combinations]
            matrix = element.matrix[elementMatrixName].reshape(-1, 1)
            matrix = [float(x) for x in matrix]
            self.matrix[systemMatrixName] += sps.coo_matrix((matrix, (rows, columns)), shape=self.matrix[systemMatrixName].shape)
        return None

    def set_vector(self, vectorName, dof, value):
        columns = [0 for x in dof]
        self.vector[vectorName] += sps.coo_matrix((value, (dof, columns)), shape=self.vector[vectorName].shape)

    def add_boundary_condition(self, boundary_condition, vectorName, matrixName):
        for bc in boundary_condition:
            if bc.type == 'load':
                for case in bc.case:
                    dof = bc.case[case][0]
                    value = bc.case[case][1]
                    self.set_vector(vectorName=vectorName, dof=dof, value=value)

            if bc.type == 'displacement':
                for case in bc.case:
                    dof = bc.case[case][0]
                    matrix = [1e20 for x in dof]
                    self.matrix[matrixName] = sps.coo_matrix((matrix, (dof, dof)), shape=self.matrix[matrixName].shape)
