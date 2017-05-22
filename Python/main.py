from NODE import NODE
from ELEMENT import ELEMENT
from SYSTEM import LinearElasticity
from BoundaryCondtion import LOAD, DISPLACEMENT

import sys
import numpy as np

grid = [[0, 0, 0], [1, 0, 0], [2, 0, 0], [3, 0, 0]]

node = []
for id in range(0, len(grid)):
    node.append(NODE())
    node[-1].id = id
    node[-1].x = grid[id][0]
    node[-1].y = grid[id][1]
    node[-1].z = grid[id][2]
    node[-1].get_dof()

element = []
for i in range(0, len(grid) - 1):
    element.append(ELEMENT())
    nodes = [node[i], node[i + 1]]
    element[-1].id = i + 1
    element[-1].attach_nodes(nodes)
    element[-1].type = 'TRUSS'
    element[-1].geometry['A'] = 10.0
    element[-1].geometry.update({'J': 1e3, 'Ix': 0.0, 'Iy': 0.0, 'Iz': 4e3})
    element[-1].material = {'E': 200e9, 'rho': 1000, 'G': 1e9}
    element[-1].initialize()

system = LinearElasticity()
system.add_elements(element)
system.add_matrix('stiffness')
system.add_vector('rhs')
system.set_matrix(systemMatrixName='stiffness', elementMatrixName='stiffness')
# print system.matrix['stiffness']
# np.savetxt('test.txt', system.matrix['stiffness'].todense(), fmt='%12.4E')

load = LOAD()
load.nodes = [0, 1]
load.magnitude = [100, 0, 0, 0, 0, 0]
load.initialize()
print load.dof
print load.value

np.savetxt('test.txt', system.matrix['stiffness'].todense(), fmt='%12.4E')