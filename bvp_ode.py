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
        self.populate_matrix_6th_order() if self.uniform else self.populate_matrix()
        self.populate_vector()
        self.apply_boundary_conditions()
        self.mtx = sparse.csr_matrix(self.mtx)
        self.mtx = self.mtx.astype(np.float64)
        self.U = dsolve.spsolve(self.mtx, np.array(self.vec), use_umfpack=True)

    def make_x_grid(self):
        self.grid_x = np.linspace(self.ode.x_min, self.ode.x_max, num=self.num_x_nodes)

    def make_y_grid(self):
        self.grid_y = np.linspace(self.ode.y_min, self.ode.y_max, num=self.num_y_nodes)


class BvpOde1D(BvpOde):

    """docstring for BvpOde1D"""

    def __init__(self, ode, bc, num_x_nodes, custom_grid=False, grid_x=False):
        self.ode = ode
        self.bc = bc
        self.num_x_nodes = num_x_nodes
        if custom_grid:
            self.grid_x = grid_x
            self.uniform = False
        else:
            self.make_x_grid()
            self.uniform = True

    def set_file_name(self, str):
        self.file_name = str

    def populate_matrix(self):
        x_vec = self.grid_x
        xm = x_vec[0:-2]
        x = x_vec[1:-1]
        xp = x_vec[2:]

        alpha = 2.0 / (xp - xm) / (x - xm)
        betta = -2.0 / (xp - x) / (x - xm)
        gamma = 2.0 / (xp - xm) / (xp - x)

        left_diag = alpha * self.ode.uxx - self.ode.ux / (xp - xm)
        main_diag = betta * self.ode.uxx + self.ode.u
        right_diag = gamma * self.ode.uxx + self.ode.ux / (xp - xm)
        main_diag = np.append([0], main_diag, 0)
        main_diag = np.append(main_diag, [0], 0)
        left_diag = np.append(left_diag, [[0, 0]])
        right_diag = np.append([0, 0], right_diag)
        self.mtx = sparse.spdiags([left_diag, main_diag, right_diag], [-1, 0, 1], self.num_x_nodes, self.num_x_nodes, format="lil")

    def populate_matrix_6th_order(self):
        h = self.grid_x[self.num_x_nodes - 1] - self.grid_x[self.num_x_nodes - 2]
        D = self.ode.uxx
        w = self.ode.ux
        k = self.ode.u
        self.mtx = sparse.lil_matrix((self.num_x_nodes, self.num_x_nodes))

        i1 = 1
        self.mtx[i1, i1 + 6] = 11 * D / (180 * h * h)
        self.mtx[i1, i1 + 5] = -90 * D / (180 * h * h) + (+2 * w) / (60 * h)
        self.mtx[i1, i1 + 4] = 324 * D / (180 * h * h) + (-15 * w) / (60 * h)
        self.mtx[i1, i1 + 3] = -670 * D / (180 * h * h) + (+50 * w) / (60 * h)
        self.mtx[i1, i1 + 2] = 855 * D / (180 * h * h) + (-100 * w) / (60 * h)
        self.mtx[i1, i1 + 1] = -486 * D / (180 * h * h) + (+150 * w) / (60 * h)
        self.mtx[i1, i1] = -70 * D / (180 * h * h) + (-77 * w) / (60 * h) + k
        self.mtx[i1, i1 - 1] = 126 * D / (180 * h * h) + (-10 * w) / (60 * h)

        i2 = 2
        self.mtx[i2, i2 + 5] = -2 * D / (180 * h * h)
        self.mtx[i2, i2 + 4] = +16 * D / (180 * h * h) + (-w) / (60 * h)
        self.mtx[i2, i2 + 3] = -54 * D / (180 * h * h) + (+8 * w) / (60 * h)
        self.mtx[i2, i2 + 2] = +85 * D / (180 * h * h) + (-30 * w) / (60 * h)
        self.mtx[i2, i2 + 1] = +130 * D / (180 * h * h) + (+80 * w) / (60 * h)
        self.mtx[i2, i2] = -378 * D / (180 * h * h) + (-35 * w) / (60 * h) + k
        self.mtx[i2, i2 - 1] = +214 * D / (180 * h * h) + (-24 * w) / (60 * h)
        self.mtx[i2, i2 - 2] = -11 * D / (180 * h * h) + (+2 * w) / (60 * h)

        im2 = self.num_x_nodes - 3
        self.mtx[im2, im2 - 5] = -2 * D / (180 * h * h)
        self.mtx[im2, im2 - 4] = +16 * D / (180 * h * h) + (+w) / (60 * h)
        self.mtx[im2, im2 - 3] = -54 * D / (180 * h * h) + (-8 * w) / (60 * h)
        self.mtx[im2, im2 - 2] = +85 * D / (180 * h * h) + (+30 * w) / (60 * h)
        self.mtx[im2, im2 - 1] = +130 * D / (180 * h * h) + (-80 * w) / (60 * h)
        self.mtx[im2, im2] = -378 * D / (180 * h * h) + (+35 * w) / (60 * h) + k
        self.mtx[im2, im2 + 1] = +214 * D / (180 * h * h) + (+24 * w) / (60 * h)
        self.mtx[im2, im2 + 2] = -11 * D / (180 * h * h) + (-2 * w) / (60 * h)

        im1 = self.num_x_nodes - 2
        self.mtx[im1, im1 - 6] = 11 * D / (180 * h * h)
        self.mtx[im1, im1 - 5] = -90 * D / (180 * h * h) + (-2 * w) / (60 * h)
        self.mtx[im1, im1 - 4] = 324 * D / (180 * h * h) + (+15 * w) / (60 * h)
        self.mtx[im1, im1 - 3] = -670 * D / (180 * h * h) + (-50 * w) / (60 * h)
        self.mtx[im1, im1 - 2] = 855 * D / (180 * h * h) + (100 * w) / (60 * h)
        self.mtx[im1, im1 - 1] = -486 * D / (180 * h * h) + (-150 * w) / (60 * h)
        self.mtx[im1, im1] = -70 * D / (180 * h * h) + (+77 * w) / (60 * h) + k
        self.mtx[im1, im1 + 1] = 126 * D / (180 * h * h) + (+10 * w) / (60 * h)

        for i in xrange(3, self.num_x_nodes - 3):
            diffusion_alpha_m3 = 2 * D / (180 * h * h) + (-w) / (60 * h)
            diffusion_alpha_m2 = -27 * D / (180 * h * h) + (+9 * w) / (60 * h)
            diffusion_alpha_m1 = 270 * D / (180 * h * h) + (-45 * w) / (60 * h)
            diffusion_alpha_0 = -490 * D / (180 * h * h) + k
            diffusion_alpha_p1 = 270 * D / (180 * h * h) + (+45 * w) / (60 * h)
            diffusion_alpha_p2 = -27 * D / (180 * h * h) + (-9 * w) / (60 * h)
            diffusion_alpha_p3 = 2 * D / (180 * h * h) + (+w) / (60 * h)

            self.mtx[i, i - 3] = diffusion_alpha_m3
            self.mtx[i, i - 2] = diffusion_alpha_m2
            self.mtx[i, i - 1] = diffusion_alpha_m1
            self.mtx[i, i] = diffusion_alpha_0
            self.mtx[i, i + 1] = diffusion_alpha_p1
            self.mtx[i, i + 2] = diffusion_alpha_p2
            self.mtx[i, i + 3] = diffusion_alpha_p3

    def populate_vector(self):
        x_vec = self.grid_x
        self.vec = [self.ode.rhs_func(x) for x in x_vec]

    def apply_boundary_conditions(self):
        D = self.ode.uxx
        w = self.ode.ux

        if self.bc.x0_is_dirichlet:
            self.mtx[0, 0] = 1
            self.vec[0] = self.bc.x0_value
            self.is_x0_bc_applied = True

        if self.bc.xn_is_dirichlet:
            self.mtx[self.num_x_nodes - 1, self.num_x_nodes - 1] = 1
            self.vec[self.num_x_nodes - 1] = self.bc.xn_value
            self.is_xn_bc_applied = True

        if self.bc.x0_is_neumann:
            h = self.grid_x[1] - self.grid_x[0]
            F = self.bc.x0_value
            self.mtx[0, 0] = -2.0 * D / (h * h)
            self.mtx[0, 1] = 2.0 * D / (h * h)
            self.vec[0] = 2 * F * (D / h - h * w)
            self.is_x0_bc_applied = True

        if self.bc.xn_is_neumann:
            h = self.grid_x[self.num_x_nodes - 1] - self.grid_x[self.num_x_nodes - 2]
            F = self.bc.xn_value
            self.mtx[self.num_x_nodes - 1, self.num_x_nodes - 2] = -2.0 * D / (h * h)
            self.mtx[self.num_x_nodes - 1, self.num_x_nodes - 1] = 2.0 * D / (h * h)
            self.vec[self.num_x_nodes - 1] = 2 * F * (D / h - h * w)
            self.is_xn_bc_applied = True

        assert self.is_x0_bc_applied and self.is_xn_bc_applied


class BvpOde2D(BvpOde1D):

    """docstring for BvpOde2D"""

    def __init__(self, ode, bc, num_x_nodes, num_y_nodes):
        BvpOde1D.__init__(self, ode, bc, num_x_nodes)
        self.num_y_nodes = num_y_nodes
        self.make_y_grid()

    def populate_matrix(self):
        pass
