
diffCoeff = 400
convCoeff = (-1.,)
# sourceCoeff = 2.
#
# We define a 1D mesh
#
# .. index:: Grid1D
#
from fipy import *
#
nx = 1000
L = 10.
mesh = Grid1D(dx=L / nx, nx=nx)
dt = 0.01

valueLeft = 0.15

var = CellVariable(name="variable", mesh=mesh)

var.constrain(valueLeft, mesh.facesLeft)
var.faceGrad.constrain([0], mesh.facesRight)


def bio():
    # print 123
    14.4 * numerix.exp(-0.25 * x)

# print var - var[0]
# raise SystemExit
# We define the convection-diffusion equation with source
eq = TransientTerm() == (DiffusionTerm(coeff=diffCoeff)
                         + PowerLawConvectionTerm(coeff=convCoeff))
                         # 2*
                         # + 14.4 * numerix.exp(-0.25 * mesh.cellCenters[0]) * (var - var[0]))
                         # + ImplicitSourceTerm(2))
#
# .. .. index:: DefaultAsymmetricSolver
#
eq.solve(var=var, dt=dt)

# axis = 0
# x = mesh.cellCenters[axis]
# AA = -sourceCoeff * x / convCoeff[axis]
# BB = 1. + sourceCoeff * L / convCoeff[axis]
# CC = 1. - numerix.exp(-convCoeff[axis] * x / diffCoeff)
# DD = 1. - numerix.exp(-convCoeff[axis] * L / diffCoeff)
# analyticalArray = AA + BB * CC / DD
viewer = Viewer(vars=var, datamin=0, datamax=0.2)

if __name__ == '__main__':
    for i in xrange(1, 100):
        eq.solve(var=var, dt=dt)
        # print var.allclose(analyticalArray, rtol=1e-4, atol=1e-4)
        viewer.plot()
        # print mesh
        # raw_input("Done")
