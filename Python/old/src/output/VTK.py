class VTK_IO():
    def __init__(self):
        self.vector_data = dict()
        self.scalar_data = dict()
        self.x = None
        self.y = None
        self.z = None
        self.points = None
        self.lines = None

    def add_points(self, points):
        self.points = points

    def add_vector(self, name, data):
        self.vector_data[name] = data

    def add_scalar(self, name, data):
        self.scalar_data[name] = data

    def add_lines(self, lines):
        self.lines = lines

    def write_to_file(self, file_path):
        self.file_path = file_path
        self.write_vtk_header('sample file', 'ASCII')
        self.write_vtk_polydata()
        self.write_vtk_lines()
        self.write_vtk_vectors()


    def write_vtk_header(self, title, data_type):
        title = title + '\n'
        data_type = data_type + '\n'
        with open(self.file_path, 'w') as vtkfile:
            vtkfile.write('# vtk DataFile Version 2.0\n')
            vtkfile.write(title)
            vtkfile.write(data_type)

    def write_vtk_polydata(self):
        with open(self.file_path, 'a') as vtkfile:
            vtkfile.write('DATASET POLYDATA\n')
            vtkfile.write('POINTS {:<6d} FLOAT\n'.format(len(self.points)))
            for point in self.points:
                vtkfile.write('{:<16.6E} {:<16.6E} {:<16.6E}\n'.format(point[0], point[1], point[2]))

    def write_vtk_lines(self):
        with open(self.file_path, 'a') as vtkfile:
            vtkfile.write('LINES {:<6d} {:<6d}\n'.format(len(self.lines), 3 * len(self.lines)))
            for line in self.lines:
                vtkfile.write('2 {:<6d} {:<6d}\n'.format(line[0], line[1]))

    def write_vtk_vectors(self):
        with open(self.file_path, 'a') as vtkfile:
            vtkfile.write('POINT_DATA {:<6d}\n'.format(len(self.points)))
            for vectors in self.vector_data:
                vtkfile.write('VECTORS {:s} float\n'.format(vectors))
                for vector in self.vector_data[vectors]:
                    vtkfile.write('{:<16.6E} {:<16.6E} {:<16.6E}\n'.format(vector[0], vector[1], vector[2]))