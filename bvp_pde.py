from bvp_ode import *

class BvpPde1D(BvpOde1D):
    """docstring for BvpPde1D"""
    def __init__(self, ode, bc, dt, tau, T, num_x_nodes, uj0, custom_grid=False, grid_x=False):
        BvpOde1D.__init__(self, ode, bc, num_x_nodes, custom_grid=False, grid_x=False)
        self.dt = dt
        self.tau = tau
        self.T = T
        self.uj0 = uj0
        self.populate_operators()
        self.populate_init_vector()
        
    def solve(self):
        self.differentiate_in_time()

    def populate_operators(self):
        I = np.identity(self.num_x_nodes)
        self.populate_matrix()
        self.populate_vector()
        self.apply_boundary_conditions()
        F =1*self.mtx
        self.operator_uj0 = (I +       self.tau * self.dt * F)
        self.operator_uj1 = (I - (1 - self.tau) * self.dt * F)

    def populate_init_vector(self):
        x_vec = self.grid_x
        self.b = np.array([self.uj0(x) for x in x_vec])
        self.Ut = np.array([self.uj0(x) for x in x_vec])

    def differentiate_in_time(self):
        for x in xrange(0,int(self.T/self.dt)):
            self.differentiate_pde_1TS()

    def differentiate_pde_1TS(self):
        self.b[0] = self.vec[0]
        self.b[self.num_x_nodes-1] = self.vec[self.num_x_nodes-1]
        temp = np.dot(np.array(self.operator_uj0), np.array(self.b))
        solution = dsolve.spsolve(sparse.csr_matrix(self.operator_uj1), np.array(temp), use_umfpack=True)
        self.b = solution
        self.Ut = np.vstack([self.Ut, solution])
        self.U = self.b
        self.b[0] = self.vec[0]
        self.b[self.num_x_nodes-1] = self.vec[self.num_x_nodes-1]