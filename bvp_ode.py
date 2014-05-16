import numpy as np
from scipy import sparse
from scipy.sparse.linalg import dsolve

class BvpOde:
    """docstring for BvpOde"""
    def __init__(self, ode, bc, num_nodes):
        self.ode = ode
        self.bc = bc
        self.num_nodes = num_nodes
        self.grid  = self.form_grid()
        self.populate_matrix()
        self.populate_vector()
        self.apply_boundary_conditions()
        self.solve()

    def form_grid(self):
        return {'x': np.linspace(self.ode.x_min, self.ode.x_max, num=self.num_nodes)}

    def set_file_name(self, str):
        self.file_name = str

    def populate_matrix(self):
        x_vec = self.grid['x']
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
        self.mtx = sparse.spdiags([left_diag, main_diag, right_diag ], [-1,0, 1], self.num_nodes, self.num_nodes, format="lil")

    def populate_vector(self):
        x_vec = self.grid['x']
        self.vec = [self.ode.rhs_func(x) for x in x_vec]

    def apply_boundary_conditions(self):
        left_bc_applied = False; 
        right_bc_applied = False;

        if self.bc.x0_bc_is_dirichlet:
            self.mtx[0,0] = 1
            left_bc_applied = True

        if self.bc.xn_bc_is_dirichlet:
            self.mtx[self.num_nodes-1,self.num_nodes-1] = 1
            right_bc_applied = True

        if self.bc.x0_bc_is_neumann:
            h = self.grid['x'][1]-self.grid['x'][0]
            self.mtx[0,0] = -1.0/h
            self.mtx[0,1] = 1/h
            left_bc_applied = True

        if self.bc.xn_bc_is_neumann:
            h = self.grid['x'][self.num_nodes-1]-self.grid['x'][self.num_nodes-2]
            self.mtx[self.num_nodes-1, self.num_nodes-2] = -1.0/h
            self.mtx[self.num_nodes-1, self.num_nodes-1] = 1/h
            right_bc_applied = True
        self.vec[0] = self.bc.x0_value
        self.vec[self.num_nodes-1] = self.bc.xn_value
        assert right_bc_applied and left_bc_applied

    def solve(self):
        self.mtx = sparse.csr_matrix(self.mtx)
        self.mtx = self.mtx.astype(np.float64)
        self.x = dsolve.spsolve(self.mtx, np.array(self.vec), use_umfpack=True)

        

        
