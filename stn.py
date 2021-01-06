import numpy as np
import igraph as ig
from scipy.integrate import odeint



def henon_map(N, a, b, x0, y0):
    xs = np.zeros(N)
    ys = np.zeros(N)
    xs[0] = x0
    ys[0] = y0

    for i in range(1,len(xs)):
        xs[i] = 1 - a * xs[i-1]**2 + ys[i-1]
        ys[i] = b * xs[i-1]

    return np.array([xs, ys])



def lorenz_step(X, t, sigma, beta, rho):
    u, v, w = X
    up = sigma*(v - u)
    vp = u*(rho-w) - v
    wp = u*v - beta*w
    return up, vp, wp



def lorenz(N=10000, dt=0.01, sigma=10, beta=8./3, rho=28, xyz0=[10,10,40]):
    t = np.arange(0, N*dt, dt) # N number of points, dt timestep
    u0, v0, w0 = xyz0  # intial conditions
    f = odeint(lorenz_step, (u0, v0, w0), t, args=(sigma, beta, rho))
    return f.T



def poincare(data_in, axis=0, value=0., positive_direction=False):
    '''
    Create poincare section of quasi-continuous timeseries in data_in.
    The plane of the section will be perpedicular to the selected axis,
    intersecting the axis at value.
    '''
    axes = list(range(len(data_in)))
    axes.remove(axis)

    d_cut  = data_in[axis]
    d_rest = data_in[axes]

    ds = np.vstack([d_cut[:-1],d_cut[1:]])
    if positive_direction:
        i  = np.where(((ds[0,:]<value)&(ds[1,:]>value)))[0]
    else:
        i  = np.where((ds[0,:]>value)&(ds[1,:]<value))[0]

    r1 = d_cut[i]
    r2 = d_cut[i+1]

    V1 = d_rest[:,i]
    V2 = d_rest[:,i+1]

    A = (V2-V1)/(r2-r1)
    B = V1-A*r1

    Vcut = A*value+B

    return Vcut


#-----------------------------------------------------------------------
def __discr__(data, b):
    a = []
    # there will be b bins between the minimum and maximum of the signal
    offset = data.min(axis=1)
    data = (data.T-offset).T
    size = data.max(axis=1) - data.min(axis=1)
    data = (data.T/size).T
    for i in range(0,len(data)):
        v = np.floor(.25+data[i]*(b-.75)*1.)
        a.append(v-v.min())
    return np.array(a,dtype=int)

#-------------------------------------------------------------------------------
def count_same(x):
  s = np.concatenate([[x.min() - 1], np.sort(x), [x.max()+1]])
  w = np.where((s[1:] - s[:-1]) > 0)[0]
  return np.vstack([s[w[1:]], w[1:] - w[:-1]])


#-----------------------------------------------------------------------
def count_links(edges, directed=True):
  """
  edges - numpy.ndarray of size (N, 2)

  Returns:
    numpy.ndarray of size (N, 3). First two columns are the different edges, third column the number of occurences
  """
  x = edges.copy()
  xmax = x.max()+1
  if not directed:
    x = np.sort(x, axis=1)
  val, count = count_same(x[:, 0]*xmax + x[:, 1])
  return np.vstack([val//xmax, val%xmax, count]).T


#-----------------------------------------------------------------------
def STN(data, b, vs_coords=True, no_dead_ends=True):
    """
    Constructs the State-Transition Network from time series data

    Parameters
    ----------
    data : numpy.ndarray
        Array with dynamics data of shape (vars, time)
    b : int
        Resolution of the spatial discretization.
        b number of bins will be constructed between the extrema of the data.
    vs_coords : bool, optional
        Whether to keep the spatial layout of the vertices.
    no_dead_ends : bool, optional
        If True then vertices with zero out-degree are removed.
    Returns
    -------
    g : igraph.Graph
        State-Transition Network.
    """
    d = __discr__(data, b)
    d = np.array(d, dtype=int)
    sh = d.max(axis=1)+1
    ts = np.ravel_multi_index(d, sh)
    esl = np.array([ts[:len(ts)-1],ts[1:]]).T
    es_ws = count_links(esl)
    #print(esl)
    g = ig.Graph((es_ws[:,:2]).tolist(),directed=True)
    if vs_coords:
        g.vs['x'], g.vs['y'] = np.unravel_index(np.arange(len(g.vs)),sh)
        g.vs['y'] = -np.array(g.vs['y'])
    g.es['weight'] = es_ws[:,2]
    zd = g.vs.select(lambda x:x.degree()==0)
    zd.delete()
    minOutDeg = lambda g:min(g.outdegree())

    # removes vertices that have zero outgoing edges
    if no_dead_ends:
        while minOutDeg(g)==0:
            nod = g.vs.select(_outdegree=0)
            g.delete_vertices(nod)
        #print 'Dead end removed!'
    return g
