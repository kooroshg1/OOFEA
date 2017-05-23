class VTK_IO():
    def __init__(self):
        self.fileName = None
        self.nodes = []
        self.connectivity = []
        self.cellType = None
        self.vector = dict()
        self.scalar = dict()
        return None

    def setFileName(self, fileName):
        self.fileName = fileName

    def setNodes(self, nodes):
        self.nodes = nodes

    def setNodeConnectivity(self, nodeConnectivity):
        self.connectivity = nodeConnectivity

    def addVector(self, vectorName):
        self.vector[vectorName] = []

    def setVector(self, vectorName, vectorValue):
        self.vector[vectorName] = vectorValue

    def addScalar(self, scalarName):
        self.scalar[scalarName] = []

    def setScalar(self, scalarName, scalarValue):
        self.scalar[scalarName] = scalarValue

    def setCellType(self, cellType):
        if cellType == 'line':
            self.cellType = 3
        elif cellType == 'triangle':
            self.cellType = 5
        elif cellType == 'quad':
            self.cellType = 9

    def write(self):
        header = """# vtk DataFile Version 2.0
{:s}
ASCII
DATASET UNSTRUCTURED_GRID\n"""
        with open(self.fileName + ".vtk", 'w') as output:
            output.write(header.format(self.fileName))
            output.write('POINTS {:<12d} FLOAT\n'.format(len(self.nodes)))
            for node in self.nodes:
                output.write('{:<16.4E} {:<16.4E} {:<16.4E}\n'.format(node[0], node[1], node[2]))
            print len(self.connectivity) + len(sum(self.connectivity, []))
            output.write('CELLS {:<16d} {:<16d}\n'.format(len(self.connectivity),
                                                          len(self.connectivity) + len(sum(self.connectivity, []))))
            for element in self.connectivity:
                output.write('{:<16d}'.format(len(element)))
                for node in element:
                    output.write(' {:<16d}'.format(node))
                output.write('\n')

            output.write('CELL_TYPES {:16d}\n'.format(len(self.connectivity)))
            for element in self.connectivity:
                output.write('{:<16d}\n'.format(self.cellType))

            output.write('POINT_DATA {:16d}\n'.format(len(self.nodes)))
            for vector in self.vector:
                output.write('VECTORS {:s} float\n'.format(vector))
                for element in self.vector[vector]:
                    output.write('{:<16.4E} {:<16.4E} {:<16.4E}\n'.format(element[0], element[1], element[2]))

