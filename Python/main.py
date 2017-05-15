import numpy as np
import scipy.sparse as sps
import scipy.sparse.linalg as spsol

class DOMAIN():
    def __init__(self, input_file_path):
        self.input_file_path = input_file_path
        self.read_mesh()

    def read_mesh(self):
        self.dimension = 6
        with open(self.input_file_path, 'r') as input_file:
            self.node = dict()
            self.element = dict()
            self.property = dict()
            self.material = dict()
            self.boundary_condition = dict()
            self.node_set = dict()
            self.element_set = dict()
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
                elif line[0:4] == 'NSET':
                    line_split = line.split(',')
                    self.node_set[int(line_split[1])] = []
                    for line_split_element in line_split[2:]:
                        self.node_set[int(line_split[1])].append(int(line_split_element))
                elif line[0:4] == 'ELSET':
                    line_split = line.split(',')
                    self.element_set[int(line_split[1])] = []
                    for line_split_element in line_split[2:]:
                        self.element_set[int(line_split[1])].append(int(line_split_element))
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
                elif line[0:12] == 'DISPLACEMENT':
                    line_split = line.split(',')
                    self.boundary_condition[line_split[0]] = [int(line_split[1]), float(line_split[2]), float(line_split[3]), float(line_split[4])]
                elif line[0:8] == 'ROTATION':
                    line_split = line.split(',')
                    self.boundary_condition[line_split[0]] = [int(line_split[1]), float(line_split[2]), float(line_split[3]), float(line_split[4])]
                elif line[0:5] == 'FORCE':
                    line_split = line.split(',')
                    self.boundary_condition[line_split[0]] = [int(line_split[1]), float(line_split[2]), float(line_split[3]), float(line_split[4])]
                elif line[0:6] == 'MOMENT':
                    line_split = line.split(',')
                    self.boundary_condition[line_split[0]] = [int(line_split[1]), float(line_split[2]), float(line_split[3]), float(line_split[4])]

class ELEMENT_BASE():
    def get_nodal_degree_of_freedom(self):
        self.nodal_degree_of_freedom = []
        for node_number in self.node:
            for dof in range(0, self.dimension):
                self.nodal_degree_of_freedom.append(self.dimension * node_number + dof)

class ELEMENT_1D_BASE(ELEMENT_BASE):
    def calculate_length(self):
        self.length = np.sqrt((self.node[1][0] - self.node[0][0]) ** 2. +
                              (self.node[1][1] - self.node[0][1]) ** 2. +
                              (self.node[1][2] - self.node[0][2]) ** 2.)

    def calculate_direction_cosines(self):
        self.alpha = (self.node[1][0] - self.node[0][0]) / self.length
        self.beta = (self.node[1][1] - self.node[0][1]) / self.length
        self.gamma = (self.node[1][2] - self.node[0][2]) / self.length

    def calculate_transformation_matrix(self):
        node[1][0] = self.node[1][0] - self.node[0][0]
        node[1][1] = self.node[1][1] - self.node[0][1]
        node[1][2] = self.node[1][2] - self.node[0][2]
        node[0] = [0., 0., 0.]

        L = np.sqrt((node[1][0] - node[0][0]) ** 2. +
                    (node[1][1] - node[0][1]) ** 2. +
                    (node[1][2] - node[0][2]) ** 2.)
        x12 = np.sqrt(node[1][0] ** 2. + node[1][1] ** 2.)
        x11 = np.sqrt(node[0][0] ** 2. + node[0][1] ** 2.)

        cosAlpha = (node[1][0] - node[0][0]) / (x12 - x11)
        sinAlpha = (node[1][1] - node[0][1]) / (x12 - x11)
        cosBeta = (x12 - x11) / L
        sinBeta = (node[1][2] - node[0][2]) / L

        LAMBDA = np.matrix([[cosAlpha * cosBeta, sinAlpha, cosAlpha * sinBeta],
                            [-sinAlpha * cosBeta, cosAlpha, -sinAlpha * sinBeta],
                            [-sinBeta, 0., cosBeta]])

        self.transformation_matrix = np.kron(np.eye(4, 4), LAMBDA)

class TRUSS(ELEMENT_1D_BASE):
    def __init__(self, node, property_list):
        self.node = node
        self.type = element_type
        self.property_list = property_list
        self.calculate_length()
        self.calculate_direction_cosines()
        self.calculate_stiffness()

    def calculate_stiffness(self):
        self.stiffness_matrix = np.matrix([[1, -1], [-1, 1]])
        self.stiffness_matrix = self.property_list['YOUNG'] * self.property_list[
            'AREA'] / self.length * self.stiffness_matrix
        R = np.matrix([[self.alpha, self.beta, self.gamma, 0., 0., 0., 0., 0., 0., 0., 0., 0.],
                       [0., 0., 0., 0., 0., 0., self.alpha, self.beta, self.gamma, 0., 0., 0.]])
        self.stiffness_matrix = R.T * self.stiffness_matrix * R

class BEAM(ELEMENT_1D_BASE):
    def __init__(self, node, properties):
        self.node = node
        self.properties = properties
        self.calculate_length()
        self.calculate_transformation_matrix()
        self.calculate_local_element_stiffness_matrix()
        self.calculate_global_element_stiffness_matrix()

    def calculate_local_element_stiffness_matrix(self):
        self.stiffness_matrix = np.zeros([12, 12])
        self.E = self.properties['YOUNG']
        self.G = self.properties['G']
        self.area = self.properties['AREA']
        self.Iy = self.properties['Iy']
        self.Iz = self.properties['Iz']
        self.J = self.properties['J']
        self.a = self.length / 2.

        self.stiffness_matrix[0, 0] = self.area * self.E / (2. * self.a)
        self.stiffness_matrix[0, 6] = -self.area * self.E / (2. * self.a)

        self.stiffness_matrix[1, 1] = 3 * self.E * self.Iz / (2. * self.a ** 3.)
        self.stiffness_matrix[1, 5] = 3 * self.E * self.Iz / (2. * self.a ** 2.)
        self.stiffness_matrix[1, 7] = -3 * self.E * self.Iz / (2. * self.a ** 3.)
        self.stiffness_matrix[1, 11] = 3 * self.E * self.Iz / (2. * self.a ** 2.)

        self.stiffness_matrix[2, 2] = 3 * self.E * self.Iy / (2. * self.a ** 3.)
        self.stiffness_matrix[2, 4] = -3 * self.E * self.Iy / (2. * self.a ** 2.)
        self.stiffness_matrix[2, 8] = -3 * self.E * self.Iy / (2. * self.a ** 3.)
        self.stiffness_matrix[2, 10] = -3 * self.E * self.Iy / (2. * self.a ** 2.)

        self.stiffness_matrix[3, 3] = self.G * self.J / (2. * self.a)
        self.stiffness_matrix[3, 9] = -self.G * self.J / (2. * self.a)

        self.stiffness_matrix[4, 4] = 2. * self.E * self.Iy / self.a
        self.stiffness_matrix[4, 8] = 3. * self.E * self.Iy / (2. * self.a ** 2.)
        self.stiffness_matrix[4, 10] = self.E * self.Iy / self.a

        self.stiffness_matrix[5, 5] = 2. * self.E * self.Iz / self.a
        self.stiffness_matrix[5, 7] = -3. * self.E * self.Iz / (2. * self.a ** 2.)
        self.stiffness_matrix[5, 11] = self.E * self.Iz / self.a

        self.stiffness_matrix[6, 6] = self.area * self.E / (2. * self.a)

        self.stiffness_matrix[7, 7] = 3 * self.E * self.Iz / (2. * self.a ** 3.)
        self.stiffness_matrix[7, 11] = -3 * self.E * self.Iz / (2. * self.a ** 2.)

        self.stiffness_matrix[8, 8] = 3 * self.E * self.Iy / (2. * self.a ** 3.)
        self.stiffness_matrix[8, 10] = 3 * self.E * self.Iy / (2. * self.a ** 2.)

        self.stiffness_matrix[9, 9] = self.G * self.J / (2. * self.a ** 3.)

        self.stiffness_matrix[10, 10] = 2 * self.E * self.Iy / self.a

        self.stiffness_matrix[11, 11] = 2 * self.E * self.Iz / self.a

        self.stiffness_matrix = self.stiffness_matrix + self.stiffness_matrix.T - np.diag(self.stiffness_matrix.diagonal())

    def calculate_global_element_stiffness_matrix(self):
        self.stiffness_matrix = self.transformation_matrix * self.stiffness_matrix

class ELEMENT(TRUSS, BEAM):
    def __init__(self, node, element_type, property_list):
        self.dimension = 6
        self.node = node
        self.type = element_type
        self.property_list = property_list
        self.calculate_stiffness()
        self.get_nodal_degree_of_freedom()

    def calculate_stiffness(self):
        if self.type == 'TRUSS':
            truss = TRUSS(self.node, self.property_list)
            self.stiffness_matrix = truss.stiffness_matrix
        elif self.type == 'BEAM':
            beam = BEAM(self.node, self.property_list)
            self.stiffness_matrix = beam.stiffness_matrix

class BOUNDARY_CONDITION():
    def apply_displacement_boundary_condition(self):
        print self.mesh.boundary_condition

class MODEL(BOUNDARY_CONDITION):
    def __init__(self, mesh):
        self.mesh = mesh
        self.stiffness_matrix = sps.coo_matrix((len(self.mesh.node) * self.mesh.dimension, len(self.mesh.node) * self.mesh.dimension))
        self.rhs = sps.coo_matrix((len(self.mesh.node) * self.mesh.dimension, 1))

    def add_matrix(self, matrix, dof):
        combinations = [(x, y) for x in dof for y in dof]
        rows = [x for (x, y) in combinations]
        columns = [y for (x, y) in combinations]
        matrix = matrix.reshape(-1, 1)
        matrix = [float(x) for x in matrix]
        self.stiffness_matrix += sps.coo_matrix((matrix, (rows, columns)), shape=self.stiffness_matrix.shape)

    def add_rhs(self, vector, dof):
        columns = [0 for y in dof]
        rows = [x for x in dof]
        vector = [float(x) for x in vector]
        self.rhs += sps.coo_matrix((vector, (rows, columns)), shape=self.rhs.shape)

    def solve(self):
        self.solution = spsol.spsolve(self.stiffness_matrix, self.rhs)





mesh = DOMAIN('sample.inp')
fea = MODEL(mesh)
print mesh.node_set
# fea.apply_displacement_boundary_condition()

output = open('matrix.txt', 'w')
for el in mesh.element:
    property_id = mesh.element[el][0]
    material_id = mesh.property[property_id][1]
    element_type = mesh.property[mesh.element[el][0]][0]

    node = dict()
    property_list = dict()

    for nd in mesh.element[el][1:]:
        node[nd] = mesh.node[nd]

    if element_type == 'TRUSS':
        property_list['YOUNG'] = mesh.material[material_id][0]
        property_list['DENSITY'] = mesh.material[material_id][1]
        property_list['AREA'] = mesh.property[property_id][2]
    elif element_type == 'BEAM':
        property_list['YOUNG'] = mesh.material[material_id][0]
        property_list['G'] = mesh.material[material_id][1]
        property_list['DENSITY'] = mesh.material[material_id][2]
        property_list['AREA'] = mesh.property[property_id][2]
        property_list['J'] = mesh.property[property_id][3]
        property_list['Ix'] = mesh.property[property_id][4]
        property_list['Iy'] = mesh.property[property_id][5]
        property_list['Iz'] = mesh.property[property_id][6]

    element = ELEMENT(node, element_type, property_list)

    fea.add_matrix(element.stiffness_matrix, element.nodal_degree_of_freedom)
    fea.add_rhs([100], [5])

    np.savetxt(output, element.stiffness_matrix, fmt='%-12.4E')
    output.write('\n')

np.savetxt(output, fea.stiffness_matrix.todense(), fmt='%-12.4E')
output.close()

# fea.solve()


class vehicle():
    def __init__(self):
        self.value = None
        self.is_second = None
        self.sold = None

    def mark_as_sold(self):
        self.sold = True

    def is_sold(self):
        print self.sold

    def get_value(self):
        print self.value

class car(vehicle):
    def __init__(self, value, is_second, sold):
        self.value = value
        self.is_second = is_second
        self.sold = sold
