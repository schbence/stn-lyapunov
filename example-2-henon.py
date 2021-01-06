import numpy as np
import igraph as ig
from stn import *
from lyapunov import lyapunov_parallel
from matplotlib import pyplot as plt



if __name__=='__main__':
    N   = 100000 # Number of sample points of the dynamics
    trans = 2000 # Number of transient samples to neglect
    dt  = 0.002 # Simulation timestep
    a_params = np.linspace(1.20, 1.23, 100) # Control parameter range for the Henon map
    lyaps = np.zeros_like(a_params)

    for i, a in enumerate(a_params):
        print(i)
        # Generating the Henon map
        dynamics = henon_map(N, a, 0.3, 1.0, 0.5)[:, trans:]
        # Constructing the STNs
        g = STN(dynamics, 20)
        # Calculating the Lyapunov measure
        lyaps[i] = lyapunov_parallel(g)

    # Plotting the Lyapunov measure in function of the 'a' control parameter
    plt.plot(a_params, lyaps)
    plt.xlabel('$a$')
    plt.ylabel('$\Lambda$')
    plt.grid()
    plt.savefig('plots/ex-2-henon-lyap.png')
