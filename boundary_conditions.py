class BoundaryConditions1D:
    """docstring for BoundaryConditions1D"""
    def __init__(self):
        self.x0_is_dirichlet = False
        self.xn_is_dirichlet = False
        self.x0_is_neumann = False
        self.xn_is_neumann = False
        self.is_x0_bc_applied = False
        self.is_xn_bc_applied = False

    def set_x0_dirichlet_bc(self, value):
        assert not self.is_x0_bc_applied, "BC was already set."
        self.x0_is_dirichlet = True
        self.x0_value = value

    def set_xn_dirichlet_bc(self, value):
        assert not self.is_xn_bc_applied, "BC was already set."
        self.xn_is_dirichlet = True
        self.xn_value = value

    def set_x0_neumann_bc(self, value):
        assert not self.is_x0_bc_applied, "BC was already set."
        self.x0_is_neumann = True
        self.x0_value = value

    def set_xn_neumann_bc(self, value):
        assert not self.is_xn_bc_applied, "BC was already set."
        self.xn_is_neumann = True
        self.xn_value = value

class BoundaryConditions2D(BoundaryConditions1D):
    """docstring for BoundaryConditions2D"""
    def __init__(self):
        BoundaryConditions1D.__init__(self)
        self.y0_is_dirichlet = False
        self.yn_is_dirichlet = False
        self.y0_is_neumann = False
        self.yn_is_neumann = False
        self.is_y0_bc_applied = False
        self.is_yn_bc_applied = False
        
    def set_y0_dirichlet_bc(self, value):
        assert not self.is_y0_bc_applied, "BC was already set."
        self.y0_is_dirichlet = True
        self.y0_value = value

    def set_yn_dirichlet_bc(self, value):
        assert not self.is_yn_bc_applied, "BC was already set."
        self.yn_is_dirichlet = True
        self.yn_value = value

    def set_y0_neumann_bc(self, value):
        assert not self.is_y0_bc_applied, "BC was already set."
        self.y0_is_neumann = True
        self.y0_value = value

    def set_yn_neumann_bc(self, value):
        assert not self.is_yn_bc_applied, "BC was already set."
        self.yn_is_neumann = True
        self.yn_value = value

        

