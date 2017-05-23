import scipy.sparse as sps
import scipy.sparse.linalg as spsolve
from BoundaryCondtion import BoundaryCondition

class LinearElasticity(BoundaryCondition):
    def __init__(self):
        self.number_of_dof = None
        self.matrix = dict()
        self.vector = dict()
        self.elements = None
        self.boundary_condition = []
        return None

    def add_elements(self, elements):
        self.elements = elements

    def add_matrix(self, matrixName):
        self.matrix[matrixName] = sps.coo_matrix((self.number_of_dof, self.number_of_dof))

    def add_vector(self, vectorName):
        self.vector[vectorName] = sps.coo_matrix((self.number_of_dof, 1))

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

    def solve(self, systemMatrixName, systemVectorName):
        self.vector['solution'] = spsolve.spsolve(self.matrix[systemMatrixName], self.vector[systemVectorName])
