class Environment(object):

    """class for storing Environment parameters"""

    def __init__(self, porosity, w_x):
        self.porosity = porosity
        self.advection_x = w_x
        self.species = {}

    def add_specie(self, specie):
        self.species[specie.name] = specie
