class TRANSLATION:
    def __init__(self, node, magnitude, name=None):
        self.node = node
        self.magnitude = magnitude
        self.dof = []
        self.name = name
        self.type = 'displacement'
        self.get_dof()

    def get_dof(self):
        dof = [0, 1, 2]
        self.dof = [i + self.node * 6 for i in dof]

class ROTATION():
    def __init__(self, node, magnitude, name=None):
        self.node = node
        self.magnitude = magnitude
        self.dof = []
        self.name = name
        self.type = 'rotation'
        self.get_dof()

    def get_dof(self):
        dof = [3, 4, 5]
        self.dof = [i + self.node * 6 for i in dof]