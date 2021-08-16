#include "../src/treerep.h"
#include "../src/graph.h"
#include <fstream>
#include <string>

void print(const Graph::wmap& weight){
	for( Graph::wmap::const_iterator it = weight.begin(); it != weight.end(); ++it){
		std::pair<int,int> edge = it->first;
		std::cout << "(" << edge.first << "," << edge.second << ") "
			<< it->second << std::endl;
	}
}

int rand_tree_test(unsigned n){
	Graph T = rand_tree(n);
	T.print();
	DistMat D = T.metric();
	std::pair<Graph,Graph::wmap> tr = treerep(D, 0.1);
	Graph G = tr.first;
	Graph::wmap W = tr.second;
	G.print();
	print(W);
	return 0;
}

int main(int argc, char* argv[]){
	rand_tree_test(64);
	return 0;
}
