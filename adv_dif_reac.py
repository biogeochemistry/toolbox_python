from fipy import *


diffCoeff = 400
convCoeff = (-5.,)
# print type(convCoeff)
nx = 1000
L = 10.
mesh = Grid1D(dx=L / nx, nx=nx)
dt = 0.001

valueLeft = 0.15

var1 = CellVariable(name="OX", mesh=mesh)
var2 = CellVariable(name="OM", mesh=mesh)

var1.constrain(valueLeft, mesh.facesLeft)
var1.faceGrad.constrain([0], mesh.facesRight)
var2.faceGrad.constrain([500], mesh.facesLeft)
var2.faceGrad.constrain([0], mesh.facesRight)

OX = var1
OM = var2

R_s = OX * OM

eq1 = TransientTerm(var=var1) == DiffusionTerm(coeff=diffCoeff,var=var1) + PowerLawConvectionTerm(coeff=convCoeff,var=var1) - R_s
print eq1
eq2 = TransientTerm(var=var2) == DiffusionTerm(coeff=5,var=var2) + PowerLawConvectionTerm(coeff=[-1],var=var2)  - R_s


# eq1.solve(var=var1, dt=dt)
# eq2.solve(var=var2, dt=dt)
eqn=eq1&eq2
viewer = Viewer(vars=var1, datamin=0, datamax=0.2)

# if __name__ == '__main__':
for i in xrange(1, 10):
    eq1.solve(var=var1, dt=dt)
    eq2.solve(var=var2, dt=dt)
    # print var2.allclose(analyticalArray, rtol=1e-4, atol=1e-4)
    viewer.plot()
    # print mesh
    # raw_input("Done")

