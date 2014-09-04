from bvp_pde import *

class CoupledPdes(object):
    """docstring for CoupledPdes"""
    def __init__(self, *arg):
        i=0
        container ={}
        for a in arg:
            container[i]=a
            i=i+1
        self.pdes = container


    def differentiate_1TS(self):
        for pde in self.pdes:
            self.pdes[pde].differentiate_1TS_pde()


    def solve(self):
        for x in xrange(0,int(self.pdes[0].T/self.pdes[0].dt)):
            self.differentiate_1TS()


        