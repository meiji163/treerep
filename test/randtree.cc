#include "../src/treerep.h"
#include <fstream>
#include <string>

int rand_tree_test(unsigned n){
	Graph T = rand_tree(n);
	T.print();
	DistMat D = T.metric();
	D.print();
	std::pair<Graph,DistMat> tr = treerep(D, 0.1);
	tr.first.print();
	tr.second.print();
	return 0;
}

int main(int argc, char* argv[]){
	rand_tree_test(6);
	//if (argc != 2){
	//	return 1;
	//}
	//DistMat M = mat_from_mtx(std::string()+argv[1]); 
	//M.print();
	

	//Graph G = graph_from_mtx(std::string()+argv[1]);
	//DistMat D = G.metric();
	//std::pair<Graph,DistMat> res= treerep(D);
	//Graph T = res.first;
	//DistMat W = res.second;
	//T.to_mtx("tree1.mtx");
	//T.print();
	//W.print();
	return 0;
}
