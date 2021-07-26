#include "treerep.h"
#include <fstream>
#include <string>

/* Example from Vincent M. Sarich "Pinniped Phylogeny," 1969.
 * He measured ``immunological distance'' between 
 * dog, bear, raccoon, weasel, sea lion, cat, and monkey.
 */
DistMat immune( {32,48,51,50,48,98,148,
					26,34,29,33,84,136,
					42,44,44,92,152,
					44,38,86,142,
					42,89,142,
					90,142,
					148},8);

int main(int argc, char* argv[]){
	if (argc != 2){
		return 1;
	}
	Graph G = load_graph(argv[1]);
	std::cout << G.size() << std::endl;
	DistMat D = G.metric();
	std::pair<Graph,DistMat> res= treerep(D);
	Graph T = res.first;
	DistMat W = res.second;
	//T.print();
	//W.print();
	return 0;
}
