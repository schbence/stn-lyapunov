import numpy as np
import igraph as ig
import matplotlib.pyplot as plt
from stn import *
from lyapunov import lyapunov_parallel_one_path
import os

def load_trials(dir, rho):
    files = os.listdir(dir)
    ext = '.gml'
    rword = 'rho='+str(rho)+'_'
    trialint = lambda s:int(s.strip(ext).split('_')[-1].split('=')[-1])
    files = list(filter(lambda s:s.endswith(ext) and (rword in s), files))
    files.sort(key=trialint)
    gs = []
    for fn in files:
        gs.append(ig.load(dir + '/' + fn))
    return gs

def apply_measure(gs, meas, attrib=None, norm=False):
    l = np.zeros(len(gs))
    for i, g in enumerate(gs):
        if norm:
            norm_weigh(g, droploops=True, global_norm=False)
        l[i] = meas(g, attribute=attrib) if attrib else meas(g)
    return l

def multi_measure(dir, rhos, trials, meas, attrib=None):
    M = np.zeros([len(rhos), len(trials)])
    rhos = list(map(lambda n:'%.3f' % n, rhos))
    for i,r in enumerate(rhos):
        print("Parameter %d, rho = %s" % (i, r))
        if len(trials) > 1:
            gs = load_trials(dir, r)
        else:
            gs = load_single_trial(dir, r, trials[0])
        M[i] = apply_measure(gs, meas, attrib=attrib)
    return M

if __name__=='__main__':
    rhos = arange(181.620, 181.750, 0.001)
    trials = np.arange(10)
    lyaps = multi_measure('./lorenz_stns/', rhos, trials, lyapunov_parallel_one_path, attrib='weight')
    lyaps_avgd = lyaps.mean(axis=1)
    plt.plot(rhos, lyaps_avgd)
