class NODE:
    def __init__(self):
        self.id = None
        self.x = None
        self.y = None
        self.z = None
        self.dof = None

    def get_dof(self):
        dof = [0, 1, 2, 3, 4, 5]
        self.dof = [i + len(dof) * self.id for i in dof]