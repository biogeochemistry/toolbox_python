from __future__ import division  # normal division
from second_order_ode import *
from boundary_conditions import *
import numpy
from bvp_ode import *
from bvp_pde import *
from coupled_pde import *
from specie_collector import *
import matplotlib.pyplot as plt
import matplotlib
from mpl_toolkits.mplot3d import Axes3D
import math
import sys
import numpy as np


def init_cond(x):
    return 0.5

x_min = 0
x_max = 20
num_x_nodes = 256
x = np.linspace(x_min, x_max, num=num_x_nodes)

a = 'ox'
b = 'om'
d_ox = 500.0
w_ox = -5
dt = 0.01
T = 1
bc_x0_type = 'Dirichlet'
bc_x0_value = 0.15
bc_xn_type = 'Neumann'
bc_xn_value = 0
init_concentrations = init_cond

species = SpecieCollector()
species.add_specie('ox', d_ox, w_ox, dt, T, bc_x0_type, bc_x0_value, bc_xn_type, bc_xn_value, init_concentrations, x_min, x_max, num_x_nodes)
species.add_specie('ox1', d_ox, w_ox, dt, T, bc_x0_type, bc_x0_value, bc_xn_type, bc_xn_value, init_concentrations, x_min, x_max, num_x_nodes)
species.add_specie('ox2', d_ox, w_ox, dt, T, bc_x0_type, bc_x0_value, bc_xn_type, bc_xn_value, init_concentrations, x_min, x_max, num_x_nodes)
print species.all
species[a] = create_single_container(d_ox, w_ox, dt, T, bc_x0_type, bc_x0_value, bc_xn_type, bc_xn_value, init_concentrations, x_min, x_max, num_x_nodes)
species[a]['pde'].solve()

print species[a]['pde'].Ut


xa = np.linspace(0, 0.5, 10)
xb = np.linspace(0.5, 1.5, 10)
xc = np.linspace(1.5, math.pi, 12)

v = np.r_[1:2:10, 11]
x1 = np.hstack((np.linspace(0, 1, 64, endpoint=False), np.linspace(1, 2, 48, endpoint=False), np.linspace(2, math.pi, 16, endpoint=True)))
print len(a)

x1 = np.logspace(0.1, 1, N, endpoint=True)
x2 = np.logspace(0.1, 1, N, endpoint=False)


y = np.zeros(128)

plt.plot(x1, y, 'o')

plt.plot(x2, y + 0.5, 'o')

plt.ylim([-0.5, 1])

plt.show()

ode1 = SecondOrderOde1D(10, -10, 0, model_prob_1, 0, 0, 15)
bc1 = BoundaryConditions1D()
bc1.set_x0_dirichlet_bc(0.15)
bc1.set_xn_neumann_bc(0)

bvp1 = BvpPde1D(ode1, bc1, 0.001, 0, 1, 120, init_cond)
# bvp1.solve_pde()

ode2 = SecondOrderOde1D(10, -10, 0, model_prob_1, 0, 0, 15)
bc2 = BoundaryConditions1D()
bc2.set_x0_dirichlet_bc(0.15)
bc2.set_xn_neumann_bc(0)

bvp2 = BvpPde1D(ode1, bc1, 0.001, 0, 1, 120, init_cond)
# bvp2.solve_pde()

cpde = CoupledPdes(bvp1, bvp2)

# print cpde.pdes
#

# a = np.array([[1, 0], [0, 1]])
# b = np.array([[4, 1], [2, 2]])

cpde.solve()
print cpde.pdes[0].Ut


x = bvp1.grid_x
X, Y = np.meshgrid(x, np.linspace(0, 10, 1001))
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(X, Y, bvp1.Ut, rstride=100, cstride=10, linewidth=0, antialiased=False)
cset = ax.contourf(X, Y, bvp1.Ut, zdir='z', offset=0)
# cset = ax.contourf(X, Y, bvp1.Ut, zdir='x', offset=-2)
# cset = ax.contourf(X, Y, bvp1.Ut, zdir='y', offset=10)

ax.set_xlabel('X')
ax.set_xlim(0, 15)
ax.set_ylabel('Y')
ax.set_ylim(0, 10)
ax.set_zlabel('Z')
ax.set_zlim(0, 0.5)

plt.show()

ode2 = SecondOrderOde1D(1, 3, -4, model_prob_2, 0, math.pi)
bc2 = BoundaryConditions1D()
bc2.set_x0_neumann_bc(-5)
bc2.set_xn_dirichlet_bc(4)
bvp2 = BvpOde1D(ode2, bc2, 64)
bvp2.solve()
# numpy.savetct(sys.stdout, bvp2.U, fmt='%.4f')

print bvp2.U

print bvp1.U
plt.plot(bvp1.grid_x, bvp1.U)
plt.show()

# print bvp.__dict__
