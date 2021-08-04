#include "../src/treerep.h"
#include "../src/graph.h"

int main(){
	for (int i=0; i<4; ++i){
		DistMat D = mat_from_mtx("../data/immune.mtx");
		D.print();
		std::pair<Graph,DistMat> tr = treerep(D,0.1);
		DistMat W  = tr.second;
		Graph T = tr.first;
		T.print();
		W.print();
		T.to_mtx("immune_tree.mtx");
	}
}
