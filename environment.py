class Environment(object):

    """class for storing Environment parameters"""

    def __init__(self, porosity, w_x, pH = 7, phreeqc_enable = False):
        self.porosity = porosity
        self.advection_x = w_x
        self.species = {}
        self.pH = pH
        self.phreeqc_enable = phreeqc_enable

    def add_species(self, species):
        self.species[species.name] = species
