from NODE import NODE
from ELEMENT import ELEMENT

import sys
import numpy as np

grid = []
for i in range(0, 2):
    grid.append([i, 0, 0])

node = []
for id in range(0, len(grid)):
    node.append(NODE())
    node[-1].id = id
    node[-1].x = grid[id][0]
    node[-1].y = grid[id][1]
    node[-1].z = grid[id][2]
    node[-1].get_dof()

element = ELEMENT()
element.attach_nodes(node)
element.type = 'TRUSS'
element.geometry = [10.0] # A
element.material = [200e9, 1000] # E, rho
element.initialize()
# print element.stiffness_matrix
np.savetxt('test.txt', element.stiffness_matrix, fmt='%12.4E')