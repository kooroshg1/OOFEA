from ELEMENT import *
from NODE import *


class PROPERTY():
    def __init__(self):
        self.type = None
        self.material_id = None
        self.gprop = None


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
        self.node = dict()
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
        self.node[int(line[1])] = NODE()
        self.node[int(line[1])].x = float(line[2])
        self.node[int(line[1])].y = float(line[2])
        self.node[int(line[1])].z = float(line[2])

    def add_element(self, line):
        self.element.append(ELEMENT(info=line))

    def add_property(self, line):
        self.property[int(line[1])] = PROPERTY()
        self.property[int(line[1])].type = line[2]
        self.property[int(line[1])].material_id = int(line[3])
        self.property[int(line[1])].gprop = [float(i) for i in line[4:]]

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
        self.assemble_mesh()

    def assemble_mesh(self):
        for element in self.element:
            print self.node
            # print self.property[element.property].type
