from __future__ import division  # normal division
from fipy import *
from matplotlib import pylab
import math


class MathModel(object):

    """all math should be here"""

    def __init__(self, environment, nx=64, Lx=10, dt=0.001, dimension='1d', ny=64, Ly=10, nz=64, Lz=10):
        self.dt = dt
        self.var = {}
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
            self.create_fipy_variable(name, specie)

    def create_fipy_variable(self, name, specie):
        self.var[name] = CellVariable(name=name, mesh=self.mesh)
        self.apply_fipy_bc(specie)

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
            self.var[specie.name].faceGrad.constrain([specie.bc_xn_value], where=self.mesh.facesRight)
        if specie.bc_x0_type == 'Flux':
            self.var[specie.name].faceGrad.constrain([specie.bc_x0_value], where=self.mesh.facesLeft)
        if specie.bc_xn_type == 'Fixed':
            self.var[specie.name].constrain(specie.bc_xn_value, where=self.mesh.facesRight)
        if specie.bc_x0_type == 'Fixed':
            self.var[specie.name].constrain(specie.bc_x0_value, where=self.mesh.facesLeft)

        if self.dimension == '2d':
            self.var[specie.name].faceGrad.constrain([0], where=self.mesh.facesTop)
            self.var[specie.name].faceGrad.constrain([0], where=self.mesh.facesBottom)

    def run(self, Time=1):
        T = 0
        veiwers = {}
        ax = {}
        plots = len(self.environment.species)
        v = 1
        # pylab.ion()
        # fig = pylab.figure()
        for name, specie in self.environment.species.iteritems():
            # ax[name] = pylab.subplot(int(math.sqrt(plots) + 1), int(math.sqrt(plots) + 1), v)
            veiwers[name] = Viewer(vars=self.var[name], datamin=0)
            v += 1
        while T < Time:
            self.eqns.solve(dt=self.dt)
            for name, var in self.var.iteritems():
                veiwers[name].plot()
            T += self.dt
