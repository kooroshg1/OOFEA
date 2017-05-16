import numpy as np

class NODE():
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z
        self.dof = None

    def get_dof(self):
        self.dof = [0, 1, 2, 3, 4, 5]
        self.dof = [i + self.node_number * 6 for i in self.dof]


class ELEMENT(NODE):
    def __init__(self):
        self.node = []
        self.node_number = -1

    def add_node(self, x, y, z):
        self.node_number += 1
        self.node.append(NODE(x, y, z))


class MESH(ELEMENT):
    def __init__(self):
        print 'Initializing mesh ..'
        self.element = []
        self.element_number = -1
        self.number_of_nodes = 0

    def add_element(self):
        self.element_number += 1
        self.element.append(ELEMENT())

    def get_number_of_nodes(self):
        for element in self.element:
            self.number_of_nodes += (element.node_number + 1)

mesh = MESH()
mesh.add_element()
mesh.add_element()

mesh.element[0].add_node(0, 0, 0)
mesh.element[0].add_node(1, 0, 0)
mesh.element[0].add_node(1, 1, 0)
mesh.element[0].add_node(0, 1, 0)

mesh.element[1].add_node(1, 0, 0)
mesh.element[1].add_node(2, 0, 0)

mesh.get_number_of_nodes()
print mesh.number_of_nodes

# for el in mesh.element:
#     print el.node_number
#     for node in el.node:
#         print node.x, node.y, node.z