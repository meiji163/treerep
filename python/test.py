import networkx as nx
import matplotlib.pyplot as plt
from trep import *
from trep_util import load_mat, to_nx_graph

if __name__ == "__main__":
    M, N = load_mat("../data/immune.mtx")
    tree, weight = treerep(M,N)

    # draw the tree with networkx
    T = to_nx_graph(weight) 
    labels = nx.get_edge_attributes(T,'weight')
    animals = {0:"dog", 1:"bear", 2:"raccoon", 3:"weasel", 4:"seal",
                5:"sea_lion", 6:"cat", 7:"monkey"}
    plt.figure(figsize=(8,8))
    nx.draw_networkx(T,
                    nodelist=[0,1,2,3,4,5,6,7],
                    labels=animals,
                    font_size=12,
                    font_weight=2)
    plt.show()
