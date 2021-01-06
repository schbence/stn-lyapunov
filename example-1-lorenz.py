import numpy as np
import igraph as ig
from stn import *
from lyapunov import lyapunov_parallel
from matplotlib import pyplot as plt

if __name__=='__main__':
    N    = 1000000 # Number of sample points of the dynamics
    trans = 2000 # Number of transient samples to neglect
    dt   = 0.002 # Simulation timestep
    xyz0 = [19, 50, 165] # Initial condition of the dynamics
    rho  = 180.7 # Control parameter of the Lorenz system
    section = 15.0 # Point of intersection of the Poincare plane with the x-axis
    b = 20 # Spatial resolution of the discretization

    # Generating the Lorenz dynamics
    dynamics = lorenz(N=N, dt=dt, sigma=10, beta=8./3, rho=rho, xyz0=xyz0)[:, trans:]
    # Taking the Poincare section
    poin = poincare(dynamics, axis=0, value=section)
    # Constructing an STN
    g = STN(poin, b)


    fig = plt.figure(figsize=(4,6))

    # Plotting the dynamics and poincare section
    plt.subplot(2,1,1)
    plt.plot(dynamics[0],dynamics[2], 'r', linewidth=0.2)
    plt.xlabel('x')
    plt.ylabel('z')
    plt.axvline(section, c='k', ls='--', alpha=0.6)
    plt.title('Lorenz system')

    plt.subplot(2,1,2)
    plt.plot(poin[0], poin[1], 'r.')
    plt.xlabel('y')
    plt.ylabel('z')
    plt.title('Poincare section')
    plt.subplots_adjust(left=0.18, bottom=0.08, top=0.94, hspace=0.3)
    plt.ion()
    plt.show()
    plt.savefig('plots/ex-1-lorenz.png')

    # Plotting the STN
    p = ig.plot(g)
    p.save('plots/ex-1-stn.png')
    p.show()
