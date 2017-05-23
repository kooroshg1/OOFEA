import sys
sys.path.insert(0, '../src')
sys.path.insert(0, '../io')

from NODE import NODE
from ELEMENT import ELEMENT
from SYSTEM import LinearElasticity
from LOAD import FORCE, MOMENT
from DISPLACEMENT import TRANSLATION, ROTATION
from VTK import VTK_IO

grid = [[720., 360., 0],
        [720., 0., 0],
        [360., 360., 0.],
        [360., 0., 0.],
        [0., 360., 0.],
        [0., 0., 0.]]

node = []
for id in range(0, len(grid)):
    node.append(NODE())
    node[-1].id = id
    node[-1].x = grid[id][0]
    node[-1].y = grid[id][1]
    node[-1].z = grid[id][2]
    node[-1].get_dof()

node_connectivity = [[4, 2],
                     [2, 0],
                     [5, 3],
                     [3, 1],
                     [3, 2],
                     [0, 1],
                     [4, 3],
                     [2, 5],
                     [2, 1],
                     [3, 0]]
element = []
element_id = 0
for connectivity in node_connectivity:
    element.append(ELEMENT())
    nodes = [node[connectivity[0]], node[connectivity[1]]]
    element[-1].id = element_id
    element[-1].attach_nodes(nodes)
    element[-1].type = 'TRUSS'
    element[-1].geometry['A'] = 1.0
    element[-1].geometry.update({'J': 1e3, 'Ix': 1.0, 'Iy': 1.0, 'Iz': 1.0})
    element[-1].material = {'E': 1e7, 'rho': 0.1, 'G': 1e9}
    element[-1].initialize()
    element_id += 1

system = LinearElasticity()
system.add_elements(element)
system.number_of_dof = len(node) * 6
system.add_matrix('stiffness')
system.add_vector('rhs')
system.set_matrix(systemMatrixName='stiffness', elementMatrixName='stiffness')

boundary_condition = []
boundary_condition.append(FORCE(node=1, magnitude=[0, -1e5, 0]))
boundary_condition.append(MOMENT(node=3, magnitude=[0, -1e5, 0]))
boundary_condition.append(TRANSLATION(node=0, magnitude=[None, None, 0]))
boundary_condition.append(TRANSLATION(node=1, magnitude=[None, None, 0]))
boundary_condition.append(TRANSLATION(node=2, magnitude=[None, None, 0]))
boundary_condition.append(TRANSLATION(node=3, magnitude=[None, None, 0]))
boundary_condition.append(TRANSLATION(node=4, magnitude=[0, 0, 0]))
boundary_condition.append(TRANSLATION(node=5, magnitude=[0, 0, 0]))
boundary_condition.append(ROTATION(node=0, magnitude=[0, 0, 0]))
boundary_condition.append(ROTATION(node=1, magnitude=[0, 0, 0]))
boundary_condition.append(ROTATION(node=2, magnitude=[0, 0, 0]))
boundary_condition.append(ROTATION(node=3, magnitude=[0, 0, 0]))
boundary_condition.append(ROTATION(node=4, magnitude=[0, 0, 0]))
boundary_condition.append(ROTATION(node=5, magnitude=[0, 0, 0]))

system.attach_boundary_condition(boundary_condition)
system.apply_boundary_condition(systemVectorName='rhs', systemMatrixName='stiffness')
system.solve(systemMatrixName='stiffness', systemVectorName='rhs')

# Output to vtk
sol = system.vector['solution']
sol = sol.reshape([12, 3])
vtk = VTK_IO()
vtk.setCellType('line')
vtk.setFileName('10bartruss')
vtk.setNodes(grid)
vtk.setNodeConnectivity(node_connectivity)
vtk.addVector('displacement')
vtk.setVector(vectorName='displacement', vectorValue=sol[0::2])
vtk.addVector('rotation')
vtk.setVector(vectorName='rotation', vectorValue=sol[1::2])
vtk.write()

