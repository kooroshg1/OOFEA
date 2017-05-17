# Add src/output to path
import sys
sys.path.insert(0, '../src/output')

import VTK as VTK

vtk = VTK.VTK_IO()
points = [[0, 0, 0], [0, 1, 0]]
displacement = [[0, 0, 0], [0, -1, 0]]
velocity = [[1, 0, 0], [2, 0, 0]]
lines = [[0, 1]]

vtk.add_points(points)
vtk.add_vector('displacement', displacement)
vtk.add_vector('velocity', velocity)
vtk.add_lines(lines)

vtk.write_to_file('myFile.vtk')