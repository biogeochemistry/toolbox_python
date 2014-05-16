import numpy as np
from scipy import sparse
from scipy.sparse.linalg import dsolve
mtx = sparse.spdiags([np.linspace(2.0, 3.0, num=50),np.linspace(2.0, 30.0, num=50), np.linspace(2.0, 3.0, num=50)], [-1,0, 1], 100, 100)
print mtx.todense()

rhs = np.array(np.linspace(2.0, 3.0, num=100))
mtx1 = mtx.astype(np.float32)
x = dsolve.spsolve(mtx1, rhs, use_umfpack=False)
print x
print "Error: ", mtx1 * x - rhs

mtx2 = mtx.astype(np.float64)
x = dsolve.spsolve(mtx2, rhs, use_umfpack=True)
print x
print "Error: ", mtx2 * x - rhs

mtx1 = mtx.astype(np.complex64)
x = dsolve.spsolve(mtx1, rhs, use_umfpack=False)
print x
print "Error: ", mtx1 * x - rhs

mtx2 = mtx.astype(np.complex128)
x = dsolve.spsolve(mtx2, rhs, use_umfpack=True)
print x
print "Error: ", mtx2 * x - rhs
