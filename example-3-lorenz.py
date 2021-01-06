import numpy as np
import igraph as ig
import matplotlib.pyplot as plt
from utility import *
from lyapunov import lyapunov_parallel



if __name__=='__main__':
    rhos = np.arange(181.620, 181.750, 0.001)
    trials = np.arange(10)
    lyaps = multi_measure('./lorenz_stns/', rhos, trials, lyapunov_parallel, attrib='weight')
    lyaps_avgd = lyaps.mean(axis=1)
    plt.plot(rhos, lyaps_avgd)
