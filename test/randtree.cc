#include "../src/treerep.h"
#include <fstream>
#include <string>

int rand_tree_test(unsigned n){
	Graph T = rand_tree(n);
	T.print();
	DistMat D = T.metric();
	D.print();
	std::pair<Graph,DistMat> tr = treerep(D, 0.1);
	DistMat W = tr.second;
	Graph G = tr.first;
	G.print();
	return 0;
}

int main(int argc, char* argv[]){
	rand_tree_test(7);
	return 0;
}
