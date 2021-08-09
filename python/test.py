import networkx as nx
import matplotlib.pyplot as plt
from trep import *
from trep_util import load_mat, to_nx_graph
import numpy as np

if __name__ == "__main__":
    M, N = load_mat("../data/immune.mtx")
    #G = load_graph("../data/bio-celegans.mtx")
    tree, weight = treerep(M,N)
    print(weight)

    # draw the tree with networkx
    T = to_nx_graph(weight) 
    nx.draw(T, with_labels=True)
    plt.show()


