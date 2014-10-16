class Specie(object):

    """Class for storing parameters of chemical species"""

    def __init__(self, name, D, bc_x0_type, bc_x0_value, reaction_term=0, aquatic=True, bc_xn_type='Flux', bc_xn_value=0):
        self.name = name
        self.DiffusionCoeff = D
        self.reaction_term = reaction_term
        self.bc_x0_type = bc_x0_type
        self.bc_x0_value = bc_x0_value
        self.is_aquatic = aquatic
        self.bc_xn_type = bc_xn_type
        self.bc_xn_value = bc_xn_value
