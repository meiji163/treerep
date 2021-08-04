from trep import *
import numpy as np

def load_mat(path):
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

if __name__ == "__main__":
    M, N = load_mat("../data/immune.mtx")
    #G = load_graph("../data/bio-celegans.mtx")
    tree, weight = treerep(M,N)
    print(weight)
    print(tree.adj_list())
