import numpy as np

class NODE():
    def __init__(self, info=None):
        self.id = int(info[1])
        self.x = float(info[2])
        self.y = float(info[3])
        self.z = float(info[4])
        self.dof = None
        self.get_dof()

    def get_dof(self):
        self.dof = [0, 1, 2, 3, 4, 5]
        self.dof = [i + self.id * 6 for i in self.dof]


class ELEMENT():
    def __init__(self, info=None):
        self.id = int(info[1])
        self.property = int(info[2])
        self.node = []
        for node in info[3:]:
            self.node.append(int(node))

class PROPERTY():
    def __init__(self, info=None):
        self.id = int(info[1])
        self.type = info[2]
        self.prop = [float(prop) for prop in info[3:]]

class MATERIAL():
    def __init__(self, info=None):
        self.id = int(info[1])
        self.prop = [float(prop) for prop in info[2:]]

class BOUNDARY_CONDITION():
    def __init__(self):
        self.type = None
        self.nset = None

class ANALYSIS_OPTION():
    def __init__(self):
        self.option = None
        self.value = None

class MESH(ELEMENT):
    def __init__(self):
        print 'Initializing mesh ..'
        self.element = []
        self.node = []
        self.boundary_condition = []
        self.analysis_option = []
        self.property = dict()
        self.material = dict()
        self.mesh_file_path = None
        self.mesh_options = {'NODE': self.add_node,
                             'ELEMENT': self.add_element,
                             'PROPERTY': self.add_property,
                             'MATERIAL': self.add_material,
                             'BC': self.add_boundary_condition,
                             'AUTOSPC': self.auto_single_point_constraint}

    def add_node(self, line):
        self.node.append(NODE(info=line))

    def add_element(self, line):
        self.element.append(ELEMENT(info=line))

    def add_property(self, line):
        self.property[int(line[1])] = PROPERTY(info=line)

    def add_material(self, line):
        self.material[int(line[1])] = MATERIAL(info=line)

    def add_boundary_condition(self, line):
        self.boundary_condition.append(BOUNDARY_CONDITION())
        self.boundary_condition[-1].type = line[1]
        self.boundary_condition[-1].nset = int(line[2])
        self.boundary_condition[-1].value = [float(i) for i in line[3:]]

    def auto_single_point_constraint(self, line):
        self.analysis_option.append(ANALYSIS_OPTION())
        self.analysis_option[-1].type = line[0]
        self.analysis_option[-1].value = line[1]


    def read_mesh(self):
        with open(self.mesh_file_path) as inputfile:
            for line in inputfile:
                if line[0] == '#':
                    continue
                line = line[:-1].split(',')
                try:
                    self.mesh_options[line[0]](line)
                except AttributeError:
                    continue
                except KeyError:
                    continue
                except NameError:
                    continue
                except TypeError:
                    continue

mesh = MESH()
mesh.mesh_file_path = 'sample.inp'
mesh.read_mesh()

# for option in mesh.analysis_option:
#     print option.type, option.value

# for bc in mesh.boundary_condition:
#     print bc.type, bc.nset, bc.value

# for node in mesh.node:
#     print node.id, node.x, node.y, node.z, node.dof

# for element in mesh.element:
#     print element.id, element.property, element.node

# for property in mesh.property:
#     print property.id, property.type, property.prop


# for node in mesh.node:
#     print node.id, node.x, node.y, node.z, node.dof
# node = globals()['NODE'](1)
