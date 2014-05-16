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

    def form_grid(self):
        return {'x': np.linspace(self.ode.x_min, self.ode.x_max, num=self.num_nodes)}

    def set_file_name(self, str):
        self.file_name = str

    def populate_matrix(self):
        x_vec = self.grid['x']
        xm = x_vec[0:-2]
        x = x_vec[1:-1]
        xp = x_vec[2:]
        print xm, x , xp


        diffusion_alpha = 2.0/(xp-xm)/(x-xm)
        diffusion_betta = -2.0/(xp-x)/(x-xm)
        diffusion_gamma = 2.0/(xp-xm)/(xp-x)
        
        main_diag = diffusion_betta*self.ode.uxx + self.ode.u
        left_diag = diffusion_alpha*self.ode.uxx - self.ode.ux/(xp-xm)
        right_diag = diffusion_gamma*self.ode.uxx + self.ode.ux/(xp-xm)
        main_diag = np.append(0, main_diag); left_diag = np.append(left_diag, 0); right_diag = np.append(0, right_diag);
        print  main_diag, left_diag, right_diag
        mtx = sparse.spdiags([left_diag,main_diag ,right_diag ], [-1,0, 1], self.num_nodes, self.num_nodes)
        print mtx.todense()

        
