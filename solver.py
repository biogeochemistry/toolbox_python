from fipy import *


class Chemical_system(object):

    """create system based on environment and species and solve it in external FiPy module"""

    def __init__(self, arg):
        self.arg = arg
