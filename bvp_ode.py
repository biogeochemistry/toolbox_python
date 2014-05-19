from grid import *
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import dsolve

class BvpOde:
    """docstring for BvpOde
    This parent class for all 1d 2d 3d bvp problems, Consist of methods which are common for all"""

    def set_file_name(self, str):
        self.file_name = str    

    def solve(self):
        self.populate_matric()
        self.populate_vector()
        self.apply_boundary_conditions()
        self.mtc = sparse.csr_matrix(self.mtc)
        self.mtc = self.mtc.astype(np.float64)
        self.U = dsolve.spsolve(self.mtc, np.array(self.vec), use_umfpack=True)

    def make_x_grid(self):
        self.grid_x  = np.linspace(self.ode.x_min, self.ode.x_max, num=self.num_x_nodes)

    def make_y_grid(self):
        self.grid_y  = np.linspace(self.ode.y_min, self.ode.y_max, num=self.num_y_nodes)

class BvpOde1D(BvpOde):
    """docstring for BvpOde1D"""
    def __init__(self, ode, bc, num_x_nodes):
        self.ode = ode
        self.bc = bc
        self.num_x_nodes = num_x_nodes
        self.make_x_grid()
        

    def set_file_name(self, str):
        self.file_name = str

    def populate_matric(self):
        x_vec = self.grid_x
        xm = x_vec[0:-2]
        x = x_vec[1:-1]
        xp = x_vec[2:]

        diffusion_alpha = 2.0/(xp-xm)/(x-xm)
        diffusion_betta = -2.0/(xp-x)/(x-xm)
        diffusion_gamma = 2.0/(xp-xm)/(xp-x)
        
        main_diag = diffusion_betta*self.ode.uxx + self.ode.u
        left_diag = diffusion_alpha*self.ode.uxx - self.ode.ux/(xp-xm)
        right_diag = diffusion_gamma*self.ode.uxx + self.ode.ux/(xp-xm)
        main_diag = np.append([0], main_diag, 0)
        main_diag = np.append(main_diag,[0], 0)
        left_diag = np.append(left_diag, [[0, 0]]); 
        right_diag = np.append([0, 0], right_diag);
        self.mtc = sparse.spdiags([left_diag, main_diag, right_diag ], [-1,0, 1], self.num_x_nodes, self.num_x_nodes, format="lil")

    def populate_vector(self):
        x_vec = self.grid_x
        self.vec = [self.ode.rhs_func(x) for x in x_vec]

    def apply_boundary_conditions(self):
        left_bc_applied = False; 
        right_bc_applied = False;

        if self.bc.x0_is_dirichlet:
            self.mtc[0,0] = 1
            left_bc_applied = True

        if self.bc.xn_is_dirichlet:
            self.mtc[self.num_x_nodes-1,self.num_x_nodes-1] = 1
            right_bc_applied = True

        if self.bc.x0_is_neumann:
            h = self.grid_x[1]-self.grid_x[0]
            self.mtc[0,0] = -1.0/h
            self.mtc[0,1] = 1.0/h
            left_bc_applied = True

        if self.bc.xn_is_neumann:
            h = self.grid_x[self.num_x_nodes-1]-self.grid_x[self.num_nodes-2]
            self.mtc[self.num_nodes-1, self.num_nodes-2] = -1.0/h
            self.mtc[self.num_nodes-1, self.num_nodes-1] = 1.0/h
            right_bc_applied = True
        self.vec[0] = self.bc.x0_value
        self.vec[self.num_x_nodes-1] = self.bc.xn_value
        assert right_bc_applied and left_bc_applied


        

        
