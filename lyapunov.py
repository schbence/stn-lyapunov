import numpy as np
import multiprocessing as mproc
from functools import partial

_var_dict = {}

def _randwalk(v0, t):
    p_adj = np.frombuffer(_var_dict['P']).reshape(_var_dict['psize'],-1)
    r_adj = np.frombuffer(_var_dict['R']).reshape(_var_dict['psize'],-1)
    res = np.zeros(t)   #array for resistivitys
    vis = np.arange(p_adj.shape[0])     #vertex indices
    np.random.seed()
    for i in range(1,t):
        v1 = np.random.choice(vis, p=p_adj[v0])
        res[i] = r_adj[v0,v1]
        v0 = v1
    return np.cumsum(res)

def _walk_diff(v0, original):
    t = _var_dict['t']
    l1 = _randwalk(v0, t)
    l2 = _randwalk(v0, t)
    return abs(l1-l2) if original else (l1-l2)**2

def _walk_one(v0):
    """
    Similar to "_walk_diff" but only works with a single path
    The 'original' has no effect
    """
    t = _var_dict['t']
    return _randwalk(v0, t)

# original
def _lslope(l, original):
    ls = l[-len(l)//2:]
    ts = np.arange(1,len(l))[-len(l)//2:]
    # return np.mean(ls/np.sqrt(ts)) if original else np.sum(ls*ts)/np.sum(ts*ts)
    return np.sum(ls*np.sqrt(ts))/np.sum(ts) if original else np.sum(ls*ts)/np.sum(ts*ts)

def pool_init(P,R,t,psize):
    _var_dict['P'] = P
    _var_dict['R'] = R
    _var_dict['t'] = t
    _var_dict['psize'] = psize

def lyapunov_parallel_one_path(g, t=1000, ens=100, original=False, attribute="weight"):
    """
    Similar to 'lyapunov_parallel' but only works with a single path and it
    measure variance


    Parameters
    ----------
    g : TYPE
        DESCRIPTION.
    t : TYPE, optional
        DESCRIPTION. The default is 1000.
    ens : TYPE, optional
        DESCRIPTION. The default is 100.
    original : TYPE, optional
        DESCRIPTION. The default is False.

    Returns
    -------
    TYPE
        DESCRIPTION.

    """
    #global pool_init
    dL = np.zeros([ens,t])

    Wadj = np.array(g.get_adjacency(attribute=attribute).data)
    Padj = (Wadj.T / Wadj.sum(axis=1)).T   # normed probabliity matrix - sum of rows = 1
    Radj = - np.log(Padj) # resistivity matrix

    psize = Padj.shape[0]
    P_raw = mproc.RawArray('d',psize*psize)
    P_np  = np.frombuffer(P_raw).reshape(psize,psize)
    np.copyto(P_np, Padj)

    R_raw = mproc.RawArray('d',psize*psize)
    R_np  = np.frombuffer(R_raw).reshape(psize,psize)
    np.copyto(R_np, Radj)

    v0s = np.random.randint(0,g.vcount(),ens)

    pool = mproc.Pool(10, initializer=pool_init, initargs=(P_raw, R_raw, t, psize))
    dL = np.array(pool.map(_walk_one, v0s))
    pool.close()
    dLm = dL.std(axis=0)

    return _lslope(dLm, True)
