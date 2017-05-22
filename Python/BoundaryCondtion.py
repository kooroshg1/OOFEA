class BASE:
    def __init__(self):
        self.type = None
        return None

    def add_case(self, caseName):
        self.case[caseName] = []

    def set_case(self, caseName, value, dof):
        self.case[caseName].append(value)
        self.case[caseName].append(dof)


class LOAD(BASE):
    def __init__(self):
        self.type = 'load'
        self.nodes = []
        self.nset = None
        self.magnitude = []
        self.value = []
        self.dof = []
        return None

    def initialize(self):
        dof = [0, 1, 2, 3, 4, 5]
        for node in self.nodes:
            self.dof.extend([i + 6 * node for i in dof])
            self.value.extend(self.magnitude)




class DISPLACEMENT(BASE):
    def __init__(self):
        self.type = 'displacement'
        self.case = dict()
        return None
