#include "../src/graph.h"

int main(){
	//Graph G = rand_tree(6);
	Graph G1 = graph_from_mtx("../data/graph1.mtx");
	G1.print();
	G1.remove_edge(1,2);
	//G1.remove_vertex(1);
	G1.print();
}
