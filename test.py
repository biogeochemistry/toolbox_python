import nose
from second_order_ode import *
from boundary_conditions import *
import numpy
from bvp_ode import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import math
import sys

def model_prob_1(x):
    return 1.0

def model_prob_2(x):
    return 34.0 * math.sin(x)

def solution_prob_1(x):
    return 0.5*x*(1-x)

def solution_prob_2(x):
    return (4*math.exp(x)+math.exp(-4*x))/(4*math.exp(math.pi)+math.exp(-4*math.pi))-5*math.sin(x)-3*math.cos(x)

def init_cond(x):
    return 0.5

def bc_constr_test():
    bc0 = BoundaryConditions1D()
    nose.tools.eq_(bc0.x0_is_dirichlet, False)
    nose.tools.eq_(bc0.xn_is_dirichlet, False)
    nose.tools.eq_(bc0.x0_is_neumann, False)
    nose.tools.eq_(bc0.xn_is_neumann, False)

def bc2D_constr_test():
    bc0 = BoundaryConditions2D()
    nose.tools.eq_(bc0.y0_is_dirichlet, False)
    nose.tools.eq_(bc0.yn_is_dirichlet, False)
    nose.tools.eq_(bc0.y0_is_neumann, False)
    nose.tools.eq_(bc0.yn_is_neumann, False)

def bc2D_dir_test():
    bc1 = BoundaryConditions2D()
    bc1.set_x0_dirichlet_bc(0.15)
    bc1.set_xn_dirichlet_bc(0)
    bc1.set_y0_dirichlet_bc(0.15)
    bc1.set_yn_dirichlet_bc(0)
    nose.tools.eq_(bc1.x0_value, 0.15)
    nose.tools.eq_(bc1.xn_value, 0)
    nose.tools.eq_(bc1.x0_is_dirichlet, True)
    nose.tools.eq_(bc1.xn_is_dirichlet, True)
    nose.tools.eq_(bc1.x0_is_neumann, False)
    nose.tools.eq_(bc1.xn_is_neumann, False)
    nose.tools.eq_(bc1.y0_is_dirichlet, True)
    nose.tools.eq_(bc1.yn_is_dirichlet, True)
    nose.tools.eq_(bc1.y0_is_neumann, False)
    nose.tools.eq_(bc1.yn_is_neumann, False)


def bc2D_dir_test():
    bc1 = BoundaryConditions2D()
    bc1.set_x0_neumann_bc(0.15)
    bc1.set_xn_neumann_bc(0)
    bc1.set_y0_neumann_bc(0.15)
    bc1.set_yn_neumann_bc(0)
    nose.tools.eq_(bc1.x0_value, 0.15)
    nose.tools.eq_(bc1.xn_value, 0)
    nose.tools.eq_(bc1.x0_is_dirichlet, False)
    nose.tools.eq_(bc1.xn_is_dirichlet, False)
    nose.tools.eq_(bc1.x0_is_neumann, True)
    nose.tools.eq_(bc1.xn_is_neumann, True)
    nose.tools.eq_(bc1.y0_is_dirichlet, False)
    nose.tools.eq_(bc1.yn_is_dirichlet, False)
    nose.tools.eq_(bc1.y0_is_neumann, True)
    nose.tools.eq_(bc1.yn_is_neumann, True)



def bc1D_dir_test():
    bc1 = BoundaryConditions1D()
    bc1.set_x0_dirichlet_bc(0.15)
    bc1.set_xn_dirichlet_bc(0)
    nose.tools.eq_(bc1.x0_value, 0.15)
    nose.tools.eq_(bc1.xn_value, 0)
    nose.tools.eq_(bc1.x0_is_dirichlet, True)
    nose.tools.eq_(bc1.xn_is_dirichlet, True)
    nose.tools.eq_(bc1.x0_is_neumann, False)
    nose.tools.eq_(bc1.xn_is_neumann, False)
    

def bc1D_neu_test():
    bc2 = BoundaryConditions1D()
    bc2.set_x0_neumann_bc(0.15)
    bc2.set_xn_neumann_bc(0)
    nose.tools.eq_(bc2.x0_value, 0.15)
    nose.tools.eq_(bc2.xn_value, 0)
    nose.tools.eq_(bc2.x0_is_neumann, True)
    nose.tools.eq_(bc2.xn_is_neumann, True)
    nose.tools.eq_(bc2.x0_is_dirichlet, False)
    nose.tools.eq_(bc2.xn_is_dirichlet, False)

def ode_test():
    ode1 = SecondOrderOde1D(10, -10, 0, model_prob_1, 0, 0 , 15) 
    nose.tools.eq_(ode1.uxx,   10)
    nose.tools.eq_(ode1.ux,   -10)
    nose.tools.eq_(ode1.u,      0)
    nose.tools.eq_(ode1.x_min,  0)
    nose.tools.eq_(ode1.x_max, 15)

def grid_test():
    ode1 = SecondOrderOde1D(10, -10, 0, model_prob_1, 0, 0 , 15)
    bc1 = BoundaryConditions1D()
    bc1.set_x0_dirichlet_bc(0.15)
    bc1.set_xn_neumann_bc(0)
    bvp1 = BvpPde1D(ode1, bc1, 0.01, 0, 10, 120, init_cond)
    nose.tools.eq_((bvp1.grid_x==np.linspace(0, 15, 120)).all(),True)
    

def bvpode_error_test():
    n = 128
    ode1 = SecondOrderOde1D(-1, 0.0, 0.0, model_prob_1, 0.0, 0 , 1)
    bc1 = BoundaryConditions1D()
    bc1.set_x0_dirichlet_bc(0)
    bc1.set_xn_dirichlet_bc(0)
    bvp1 = BvpOde1D(ode1, bc1, n)
    bvp1.solve()
    sum = 0
    for i in xrange(0,128):
        x = bvp1.grid_x[i]
        sum += math.pow((bvp1.U[i] - solution_prob_1(x)),2)
    sum/=n-1
    print '\nFor problem 1 the variance is: ', sum
    nose.tools.eq_((sum<math.pow(10,-16)),True)

def bvpode2_error_test():
    n = 128
    ode1 = SecondOrderOde1D(1, 3.0, -4.0, model_prob_2, 0.0, 0 , math.pi)
    bc1 = BoundaryConditions1D()
    bc1.set_x0_neumann_bc(-5.0)
    bc1.set_xn_dirichlet_bc(4.0)
    bvp1 = BvpOde1D(ode1, bc1, n)
    bvp1.solve()
    sum = 0    
    for i in xrange(0,128):
        x = bvp1.grid_x[i]
        sum += math.pow((bvp1.U[i] - solution_prob_2(x)),2)
    sum/=(n-1)
    print '\nFor problem 2 the variance is: ', sum
    nose.tools.eq_((sum<math.pow(10,-4)),True)
    
    