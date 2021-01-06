import numpy as np
import igraph as ig
from stn import *
from lyapunov import lyapunov_parallel

if __name__=='__main__':
    N   = 1000000
    trans = 2000
    dt  = 0.002
    xyz0 = [19, 50, 165]
    a_params = np.linspace(1.20, 1.23, 100)
    lyaps = np.zeros_like(a_params)

    for i, a in enumerate(a_params):
        print(i)
        lorenz_data = lorenz(N=N, dt=dt, sigma=10, beta=8./3, rho=rho, xyz0=xyz0)[:, trans:]
        poin = poincare(dynamics, cut_by=0, cut_at=15.0, positive_direction=False)
        g = STN(dynamics, 20, 1)
        lyaps[i] = lyapunov_parallel(g)
