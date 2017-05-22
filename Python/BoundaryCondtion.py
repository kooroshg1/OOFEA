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
        self.case = dict()
        return None


class DISPLACEMENT(BASE):
    def __init__(self):
        self.type = 'displacement'
        self.case = dict()
        return None
