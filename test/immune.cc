#include "../src/treerep.h"
#include "../src/graph.h"

void print(const Graph::wmap& weight){
	for( Graph::wmap::const_iterator it = weight.begin(); it != weight.end(); ++it){
		std::pair<int,int> edge = it->first;
		std::cout << "(" << edge.first << "," << edge.second << ") "
			<< it->second << std::endl;
	}
}


int main(){
	DistMat D = mat_from_mtx("../data/immune.mtx");
	for (int i=0; i<2; ++i){
		std::pair<Graph,Graph::wmap> tr = treerep(D,0.1);
		Graph T = tr.first;
		Graph::wmap W = tr.second;
		T.print();
		print(W);
		//double ad= avg_distortion(D, W);
		//std::cout << "Avg distortion: " << ad << std::endl;
	}
}
