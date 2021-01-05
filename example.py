import numpy as np
import igraph as ig
from stn import *
from lyapunov import lyapunov_parallel_one_path

if __name__=='__main__':
    N   = 2500000
    trans = 20000
    dt  = 0.002
    xyz0 = [19, 50, 165]
    rhos = np.linspace(181.65, 181.75, 50)
    lyaps = np.zeros_like(rhos)

    for i, rho in enumerate(rhos):
        print(i)
        lorenz_data = lorenz(N=N, dt=dt, sigma=10, beta=8./3, rho=rho, xyz0=xyz0)[:, trans:]
        poin = poincare(lorenz_data, cut_by=0, cut_at=15.0, positive_direction=False)
        g = STN(poin, 20, 1, vs_coords=True, discr_mode='minmax', no_dead_ends=True)
        lyaps[i] = lyapunov_parallel_one_path(g)
