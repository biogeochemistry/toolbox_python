from second_order_ode import *
from boundary_conditions import *
from bvp_ode import *

def model_prob_1(x):
    return 1.0

ode = SecondOrderOde(-1, 0, 0, model_prob_1, 0 , 1)
bc = BoundaryConditions()
bc.set_x0_dirichlet_bc(0)
bc.set_xn_dirichlet_bc(0)

bvp = BvpOde(ode, bc, 4)
bvp.set_file_name('att.dat')

# print bvp.__dict__


