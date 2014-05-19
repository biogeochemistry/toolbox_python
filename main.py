from second_order_ode import *
from boundary_conditions import *
import numpy
from bvp_ode import *
import matplotlib.pyplot as plt
import math
import sys



def model_prob_2d(x,y):
    return -4*(1-x**2-y**2)*math.exp(-(x**2+y**2))

def bc_at_x0(y):
    return math.exp(-y**2)

def bc_at_xn(y):
    return math.exp(-1-y**2)

# def bc

# uxx, ux, uyy, uy, u, rhs_func, x_min, x_max, y_min, y_max)
ode2d = SecondOrderOde2D(1, 0, 1, 0, 0, model_prob_2d, 0, 1, 0, 2)
bc2d = BoundaryConditions2D()

















# def model_prob_1(x):
#     return 1.0
# def model_prob_2(x):
#     return 34.0 * math.sin(x)

# def solution_prob_1(x):
#     return 1.0/2*x*(1-x)

# def solution_prob_2(x):
#     return (4*math.exp(x)+math.exp(-4*x))/(4*math.exp(math.pi)+math.exp(-4*math.pi))-5*math.sin(x)-3*math.cos(x)

# ode1 = SecondOrderOde1D(-1, 0, 0, model_prob_1, 0 , 1)
# bc1 = BoundaryConditions1D()
# bc1.set_x0_dirichlet_bc(0)
# bc1.set_xn_dirichlet_bc(0)
# bvp1 = BvpOde1D(ode1, bc1, 100)
# bvp1.solve()


# ode2 = SecondOrderOde1D(1, 3, -4, model_prob_2, 0, math.pi)
# bc2 = BoundaryConditions1D()
# bc2.set_x0_neumann_bc(-5)
# bc2.set_xn_dirichlet_bc(4)
# bvp2 = BvpOde1D(ode2, bc2, 64)
# bvp2.solve()
# # numpy.savetct(sys.stdout, bvp2.U, fmt='%.4f')


# print bvp2.U

# plt.plot(bvp1.grid_x,bvp1.U)
# plt.plot(bvp1.grid_x, [solution_prob_1(x) for x in bvp1.grid_x], 'r+')
# plt.show()

# print bvp.__dict__


