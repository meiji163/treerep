#include "treerep.h"
#include <fstream>
#include <string>

int main(int argc, char* argv[]){
	if (argc != 2){
		return 1;
	}
	DistMat M = mat_from_mtx(std::string()+argv[1]); 
	M.print();
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
