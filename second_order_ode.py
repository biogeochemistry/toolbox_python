class SecondOrderOde1D:
    """SecondOrderOde1D class"""
    def __init__(self, uxx, ux, u, rhs_func, sources_sinks_func, x_min, x_max):
        self.uxx = uxx
        self.ux = ux
        self.u = u
        self.rhs_func = rhs_func
        self.sources_sinks_func = sources_sinks_func
        self.x_min = x_min
        self.x_max = x_max

class SecondOrderOde2D(SecondOrderOde1D):
    """docstring for SecondOrderOde1D"""
    def __init__(self, uxx, ux, uyy, uy, u, rhs_func, sources_sinks_func, x_min, x_max, y_min, y_max):
        SecondOrderOde1D.__init__(self, uxx, ux, u, rhs_func, sources_sinks_func, x_min, x_max)
        self.uyy = uyy
        self.uy = uy
        self.y_min = y_min
        self.y_max = y_max
        