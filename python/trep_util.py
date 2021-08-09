import numpy as np
def load_mat(path):
    '''load distance matrix from a mtx file'''
    with open(path,'r') as f:
        L = f.readline()
        while L[0] == '%':
            L = f.readline()
        _, N, S = [int(i) for i in L.split(" ")]
        assert N*(N-1)//2 == S
        dist = np.zeros(S)
        for _ in range(S):
            L = f.readline()
            vals = L.split(" ")
            i, j = [int(k) for k in vals[:2]]
            val = float(vals[-1])
            i -= 1
            j -= 1
            if i > j: i,j = j,i
            dist[N*i + j - ((i+2)*(i+1))//2]=val
    return dist, N

def to_nx_graph(weights):
    '''convert dict of edges and weights to networkx Graph'''
    import networkx as nx
    G = nx.Graph()
    for e in weights:
        if weights[e] >= 0:
            G.add_edge(e[0],e[1], weight=weights[e])
    return G

def max_min(A, B):
    '''max-min product of two square matrices
    params:
        A, B: NxN numpy arrays '''
    return np.max(np.minimum(A[:, :, None], B[None, :, :]), axis=1)

def gromov_prod(x,y,z,dist):
    return 0.5*(dist[x,z]+dist[y,z]-dist[x,y])

def mat_gromov_prod(dists, base):
    '''Gromov products of N-point metric space relative to base point
    params:
        dists: NxN numpy array of pairwise distances
        base: index of the basepoint in 0...n-1 '''
    row = dists[base, :][np.newaxis, :]
    col = dists[:, base][:, np.newaxis]
    return 0.5*(row+col-dists)

def delta_rel(dists, base):
    G = mat_gromov_prod(dists, base)
    delta = np.max(max_min(G,G)-G)
    diam = np.max(dists)
    return delta/diam
