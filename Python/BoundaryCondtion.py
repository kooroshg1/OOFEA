import scipy.sparse as sps

class BoundaryCondition:
    def attach_boundary_condition(self, boundary_condition):
        try:
            self.boundary_condition.extend(boundary_condition)
        except TypeError:
            self.boundary_condition.append(boundary_condition)

    def apply_boundary_condition(self, systemVectorName, systemMatrixName):
        for bc in self.boundary_condition:
            if (bc.type == 'force') or (bc.type == 'moment'):
                self.set_vector(vectorName=systemVectorName, dof=bc.dof, value=bc.magnitude)
            elif (bc.type == 'displacement') or (bc.type == 'rotation'):
                matrix = [1e20 for x in bc.magnitude if x != None]
                el = 0
                dof = []
                magnitude = []
                for i in bc.dof:
                    if bc.magnitude[el] != None:
                        dof.append(bc.dof[el])
                        magnitude.append(bc.magnitude[el])
                    el += 1
                self.matrix[systemMatrixName] += sps.coo_matrix((matrix, (dof, dof)), shape=self.matrix[systemMatrixName].shape)
                self.set_vector(vectorName=systemVectorName, dof=dof, value=magnitude)

