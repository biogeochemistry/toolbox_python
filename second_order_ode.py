class SecondOrderOde:
    """SecondOrderOde class"""
    def __init__(self, uxx, ux, u, rhs_func, x_min, x_max):
        self.uxx = uxx
        self.ux = ux
        self.u = u
        self.rhs_func = rhs_func
        self.x_min = x_min
        self.x_max = x_max
