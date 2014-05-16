from second_order_ode import *
from boundary_conditions import *
from bvp_ode import *
import matplotlib.pyplot as plt
import math

def model_prob_1(x):
    return 1.0
def model_prob_2(x):
    return 34* math.sin(x)
# ode = SecondOrderOde(-1, 0, 0, model_prob_1, 0 , 1)
# bc = BoundaryConditions()
# bc.set_x0_dirichlet_bc(0)
# bc.set_xn_dirichlet_bc(0)
# bvp = BvpOde(ode, bc, 1000)
# bvp.set_file_name('att.dat')

for i in xrange(1,1000):
    ode2 = SecondOrderOde(1, 3, -4, model_prob_2, 0, math.pi)
    bc2 = BoundaryConditions()
    bc2.set_x0_neumann_bc(-5)
    bc2.set_xn_dirichlet_bc(4)
    bvp2 = BvpOde(ode2,bc2, 1000)

# print bvp2.x

# plt.plot(bvp2.grid['x'],bvp2.x)
# plt.show()

# print bvp.__dict__


