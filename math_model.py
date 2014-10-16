from __future__ import division  # normal division
from fipy import *



class MathModel(object):

    """all math should be here"""

    def __init__(self, environment, nx=256, Lx=10, dt=0.001):
        self.nx = nx
        self.Lx = Lx
        self.mesh1d = Grid1D(dx=Lx / nx, nx=nx)
        self.dt = dt
        self.var = {}
        self.eqns = 0
        self.environment = environment
        self.create_model()

    def create_model(self):
        self.create_variables()
        self.create_pdes()

    def create_variables(self):
        for name, specie in self.environment.species.iteritems():
            self.create_fipy_variable(name, specie)

    def create_fipy_variable(self, name, specie):
        self.var[name] = CellVariable(name=name, mesh=self.mesh1d)
        self.apply_fipy_bc(specie)

    def create_pdes(self):
        for name, specie in self.environment.species.iteritems():
            self.making_fipy_pde(name, specie)

    def making_fipy_pde(self, name, specie):
        # add equations
        self.eqns &= TransientTerm(var=self.var[name]) == (DiffusionTerm(coeff=[specie.DiffusionCoeff], var=self.var[name]) + PowerLawConvectionTerm([self.environment.advection_x], var=self.var[name]))

    def apply_fipy_bc(self, specie):
        if specie.bc_xn_type == 'Flux':
            self.var[specie.name].faceGrad.constrain([specie.bc_xn_value], where=self.mesh1d.facesRight)
        if specie.bc_x0_type == 'Flux':
            self.var[specie.name].faceGrad.constrain([specie.bc_x0_value], where=self.mesh1d.facesLeft)
        if specie.bc_xn_type == 'Fixed':
            self.var[specie.name].constrain(specie.bc_xn_value, where=self.mesh1d.facesRight)
        if specie.bc_x0_type == 'Fixed':
            self.var[specie.name].constrain(specie.bc_x0_value, where=self.mesh1d.facesLeft)


    def run(self, Time=1):
        T=0
        veiwers = {}
        for name, specie in self.environment.species.iteritems():
            veiwers[name] = Viewer(vars=self.var[name], datamin=0)
        while T < Time:
            self.eqns.solve(dt=self.dt)
            for name, var in self.var.iteritems():
                veiwers[name].plot()
            T+=self.dt
