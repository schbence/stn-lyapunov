import numpy as np
import igraph as ig
from stn import *
from lyapunov import *

N   = 250000
trans = 20000
dt  = 0.002
rho = 180.10

xyz0 = [19, 50, 165]
lorenz_data = lorenz(N=N, dt=dt, sigma=10, beta=8./3, rho=rho, xyz0=xyz0)[:, trans:]

poin = poincare(lorenz_data, cut_by=0, cut_at=15.0, positive_direction=False)
g = STN(poin, 20, 1, vs_coords=True, discr_mode='minmax', no_dead_ends=True)
l = lyapunov_parallel_one_path(g)
print(l)
