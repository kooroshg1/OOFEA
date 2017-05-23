from VTK import VTK_IO


grid = [[720., 360., 0],
        [720., 0., 0],
        [360., 360., 0.],
        [360., 0., 0.],
        [0., 360., 0.],
        [0., 0., 0.]]

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

dispalcement = [[-0.1, -1, 0],
                [-0.2, -3, 0],
                [0., -2, 0],
                [0, -3, 0],
                [0, 0, 0],
                [0, 0, 0]]

vtk = VTK_IO()
vtk.setFileName('10bartruss')
vtk.setCellType('line')
vtk.setNodes(grid)
vtk.setNodeConnectivity(node_connectivity)
vtk.addVector('displacement')
vtk.setVector('displacement', dispalcement)
vtk.write()