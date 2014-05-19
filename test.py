from bvp_ode import *
from second_order_ode import *
from boundary_conditions import *
import math
import numpy
import matplotlib.pyplot as plt

def solution_prob_1(x):
    return 1.0/2*x*(1-x)

def solution_prob_2(x):
    return (4.0*math.exp(x)+math.exp(-4.0*x))/(4.0*math.exp(math.pi)+math.exp(-4.0*math.pi))-5.0*math.sin(x)-3.0*math.cos(x)

def model_prob_1(x):
    return 1.0

def model_prob_2(x):
    return 34* math.sin(x)

# test: before all
# 
ode1 = SecondOrderOde1D(-1, 0, 0, model_prob_1, 0 , 1)
bc1 = BoundaryConditions1D()
bc1.set_x0_dirichlet_bc(0)
bc1.set_xn_dirichlet_bc(0)
bvp1 = BvpOde1D(ode1, bc1, 100)
bvp1.solve()

x_min2 = 0
x_max2 = math.pi
nodes2 = 10
uxx2 = 1
ux2 = 3
u2 = -4
bcx02 = -5
bcxn2 = 4

ode2 = SecondOrderOde1D(uxx2, ux2, u2, model_prob_2, x_min2, x_max2)
bc2 = BoundaryConditions1D()
bc2.set_x0_neumann_bc(bcx02)
bc2.set_xn_dirichlet_bc(bcxn2)
bvp2 = BvpOde1D(ode2, bc2, nodes2)
bvp2.solve()



class Test1D:

    def test_second_order_ode_1D(self):
        assert ode2.uxx == uxx2
        assert ode2.ux == ux2
        assert ode2.u == u2

    def test_boundary_conditions_init_1D(self):
        bc0 = BoundaryConditions1D()
        assert not bc0.x0_is_dirichlet, 'created BC should be False'
        assert not bc0.xn_is_dirichlet, 'created BC should be False'
        assert not bc0.x0_is_neumann, 'created BC should be False'
        assert not bc0.xn_is_neumann, 'created BC should be False'

    def test_boundary_conditions_1D(self):
        assert not bc2.x0_is_dirichlet, 'not assigned BC'
        assert bc2.xn_is_dirichlet, 'not assigned BC'
        assert bc2.x0_is_neumann, 'not assigned BC'
        assert not bc2.xn_is_neumann, 'not assigned BC'
        assert bc2.x0_value == bcx02, 'not assigned BC'
        assert bc2.xn_value == bcxn2, 'not assigned BC'

    def test_grid_formation_1d(self):
        numpy.alltrue(bvp2.grid_x == np.linspace(x_min2, x_max2, num=nodes2))

    def test_error_2(self):
        exact_sol_2 = [solution_prob_2(x) for x in bvp2.grid_x]
        assert math.sqrt(sum((bvp2.U - exact_sol_2)**2)) < 1
        exact_sol_1 = [solution_prob_1(x) for x in bvp1.grid_x]
        assert math.sqrt(sum((bvp1.U - exact_sol_1)**2)) < 1e-5

# plt.plot(bvp1.grid_x,bvp1.U)
# plt.plot(bvp1.grid_x, [solution_prob_1(x) for x in bvp1.grid_x], 'r+')
# plt.show()

class Test2D:

    def test_second_order_ode_2D(self):
        ode2d = SecondOrderOde2D(1, 2, 3, 4, 5, 6, 7, 8, 9, 0)
        assert ode2d.uxx == 1
        assert ode2d.ux == 2
        assert ode2d.uyy == 3
        assert ode2d.uy == 4
        assert ode2d.u == 5
        assert ode2d.rhs_func == 6
        assert ode2d.x_min == 7
        assert ode2d.x_max == 8
        assert ode2d.y_min == 9
        assert ode2d.y_max == 0

    def test_boundary_conditions_2D(self):
        bc0 = BoundaryConditions2D()
        assert not bc0.x0_is_dirichlet, 'created BC should be False'
        assert not bc0.xn_is_dirichlet, 'created BC should be False'
        assert not bc0.x0_is_neumann, 'created BC should be False'
        assert not bc0.xn_is_neumann, 'created BC should be False'
        assert not bc0.y0_is_dirichlet, 'created BC should be False'
        assert not bc0.yn_is_dirichlet, 'created BC should be False'
        assert not bc0.y0_is_neumann, 'created BC should be False'
        assert not bc0.yn_is_neumann, 'created BC should be False'




