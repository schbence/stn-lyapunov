import numpy as np
import igraph as ig
from stn import *
from lyapunov import lyapunov_parallel
from matplotlib import pyplot as plt



if __name__=='__main__':
    N   = 100000
    trans = 2000
    dt  = 0.002
    xyz0 = [19, 50, 165]
    a_params = np.linspace(1.20, 1.23, 100)
    lyaps = np.zeros_like(a_params)

    for i, a in enumerate(a_params):
        print(i)
        dynamics = henon_map(N, a, 0.3, 1.0, 0.5)[:, trans:]
        g = STN(dynamics, 20)
        lyaps[i] = lyapunov_parallel(g)

    plt.plot(a_params, lyaps)
    plt.xlabel('$a$')
    plt.ylabel('$\Lambda$')
    plt.grid()
    plt.savefig('plots/ex-2-henon-lyap.png')
