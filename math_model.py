from __future__ import division  # normal division
from fipy import *
from matplotlib import pylab
import math
import species
# import pdb
from IPython.core.debugger import Tracer
# from fipy.tools import parallel
import numpy as np


class MathModel(object):

    """all math should be here"""

    def __init__(self, environment, nx=64, Lx=10, dt=0.1, dimension='1d', ny=64, Ly=10, nz=64, Lz=10):
        self.dt = dt
        self.var = {}
        self.dcdt = {}
        self.rate_expr = {}
        self.eqns = 0
        self.dimension = dimension
        if dimension == '1d':
            self.mesh = Grid1D(dx=Lx / nx, nx=nx)
        elif dimension == '2d':
            self.mesh = Grid2D(dx=Lx / nx, dy=Ly / ny, nx=nx, ny=ny)
        elif dimension == '3d':
            self.mesh = Grid3D(nx=50, ny=100, nz=10, dx=Lx / nx, dy=Ly / ny, dz=Lz / nz)
        self.environment = environment
        self.create_model()

    def create_model(self):
        self.create_variables()
        self.create_pdes()

    def create_variables(self):
        for name, specie in self.environment.species.iteritems():
            self.create_fipy_variables(name, specie)
            self.populate_other_variables(name, specie)

    def create_fipy_variables(self, name, specie):
        self.var[name] = CellVariable(name=name, mesh=self.mesh)
        self.apply_fipy_bc(specie)

    def populate_other_variables(self, name, specie):
        self.dcdt[name] = np.zeros(np.shape(self.var[name].value))
        self.rate_expr[name] = specie.reaction_term

    def create_pdes(self):
        for name, specie in self.environment.species.iteritems():
            self.making_fipy_pde(name, specie)

    def making_fipy_pde(self, name, specie):
        D = specie.DiffusionCoeff
        w = self.environment.advection_x
        if self.dimension == '1d':
            # & - add equation
            self.eqns &= TransientTerm(var=self.var[name]) == (DiffusionTerm(coeff=D, var=self.var[name]) + PowerLawConvectionTerm(coeff=(w,), var=self.var[name]))
        elif self.dimension == '2d':
            # & - add equation
            self.eqns &= TransientTerm(var=self.var[name]) == (DiffusionTerm(coeff=D, var=self.var[name]) + PowerLawConvectionTerm(coeff=[w, 0], var=self.var[name]))

    def apply_fipy_bc(self, specie):
        if specie.bc_xn_type == 'Flux':
            self.var[specie.name].faceGrad.constrain([-specie.bc_xn_value], where=self.mesh.facesRight)
        if specie.bc_x0_type == 'Flux':
            self.var[specie.name].faceGrad.constrain([-specie.bc_x0_value], where=self.mesh.facesLeft)
        if specie.bc_xn_type == 'Fixed':
            self.var[specie.name].constrain(specie.bc_xn_value, where=self.mesh.facesRight)
        if specie.bc_x0_type == 'Fixed':
            self.var[specie.name].constrain(specie.bc_x0_value, where=self.mesh.facesLeft)

        if self.dimension == '2d':
            self.var[specie.name].faceGrad.constrain([0], where=self.mesh.facesTop)
            self.var[specie.name].faceGrad.constrain([0], where=self.mesh.facesBottom)

    def run(self, Time=1, graphs=False):
        T = 0
        if graphs:
            veiwers = {}
            for name, specie in self.environment.species.iteritems():
                veiwers[name] = Viewer(vars=self.var[name], datamin=0)
        while T < Time:
            self.eqns.solve(dt=self.dt)
            # for name, specie in self.environment.species.iteritems():
                # specie.conc = self.var[specie.name].value
            # self.rates()
            if graphs:
                for name, var in self.var.iteritems():
                    veiwers[name].plot()
            T += self.dt

    def rates(self):
        C0 = np.array([self.var[name].value for name, specie in self.environment.species.iteritems()])
        self.dcdt = np.zeros(np.shape(np.array(C0)))

def butcher5(C0, rates, ts):
    """Numerical integration of ODE using Butcher's Fifth-Order Runge-Kutta method."""
    dt = 1.0 / ts
    for i in xrange(0, ts):
        k_1 = dt * rates(C0)
        k_2 = dt * rates(C0 + 1 / 4 * k_1)
        k_3 = dt * rates(C0 + 1 / 8 * k_1 + 1 / 8 * k_2)
        k_4 = dt * rates(C0 - 1 / 2 * k_2 + k_3)
        k_5 = dt * rates(C0 + 3 / 16 * k_1 + 9 / 16 * k_4)
        k_6 = dt * rates(C0 - 3 / 7 * k_1 + 2 / 7 * k_2 + 12 / 7 * k_3 - 12 / 7 * k_4 + 8 / 7 * k_5)
        C_new = C0 + (7 * k_1 + 32 * k_3 + 12 * k_4 + 32 * k_5 + 7 * k_6) / 90
        C0 = C_new
    return C_new


def rk4(C0, rates, ts):
    """Numerical integration of ODE using Fourth-Order Runge-Kutta method."""
    dt = 1.0 / ts
    for i in xrange(0, ts):
        k_1 = dt * rates(C0)
        k_2 = dt * rates(C0 + 0.5 * k_1)
        k_3 = dt * rates(C0 + 0.5 * k_2)
        k_4 = dt * rates(C0 + k_3)
        C_new = C0 + (k_1 + 2 * k_2 + 2. * k_3 + k_4) / 6
        C0 = C_new
    return C_new