from __future__ import division # normal division
from second_order_ode import *
from boundary_conditions import *
import numpy
from bvp_ode import *
import matplotlib.pyplot as plt
import matplotlib
from mpl_toolkits.mplot3d import Axes3D
import math
import sys
import numpy as np


def model_prob_1(x):
    return 0

def model_prob_2(x):
    return 34.0 * math.sin(x)

def solution_prob_1(x):
    return 1.0/2*x*(1-x)

def solution_prob_2(x):
    return (4*math.exp(x)+math.exp(-4*x))/(4*math.exp(math.pi)+math.exp(-4*math.pi))-5*math.sin(x)-3*math.cos(x)

def init_cond(x):
    return 0.5




# ode1 = SecondOrderOde1D(10, -10, 0, model_prob_1, 0, 0 , 15)
# bc1 = BoundaryConditions1D()
# bc1.set_x0_dirichlet_bc(0.15)
# bc1.set_xn_neumann_bc(0)

# bvp1 = BvpPde1D(ode1, bc1, 0.01, 0, 10, 120, init_cond)
# bvp1.solve_pde()

# x = bvp1.grid_x
# X,Y = np.meshgrid(x, np.linspace(0,10,1001))
# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# ax.plot_surface(X, Y, bvp1.Ut,rstride=100, cstride=10, linewidth=0, antialiased=False)
# cset = ax.contourf(X, Y, bvp1.Ut, zdir='z', offset=0)
# # cset = ax.contourf(X, Y, bvp1.Ut, zdir='x', offset=-2)
# # cset = ax.contourf(X, Y, bvp1.Ut, zdir='y', offset=10)

# ax.set_xlabel('X')
# ax.set_xlim(0, 15)
# ax.set_ylabel('Y')
# ax.set_ylim(0, 10)
# ax.set_zlabel('Z')
# ax.set_zlim(0, 0.5)

# plt.show()

# ode2 = SecondOrderOde1D(1, 3, -4, model_prob_2, 0, math.pi)
# bc2 = BoundaryConditions1D()
# bc2.set_x0_neumann_bc(-5)
# bc2.set_xn_dirichlet_bc(4)
# bvp2 = BvpOde1D(ode2, bc2, 64)
# bvp2.solve()
# # numpy.savetct(sys.stdout, bvp2.U, fmt='%.4f')

# print bvp2.U

# print bvp1.U
# plt.plot(bvp1.grid_x, bvp1.U)
# plt.show()

# print bvp.__dict__




