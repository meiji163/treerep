#include <iostream>
#include "treerep.hpp"

int main(){

	std::cout << "Test 1" << std::endl;
	/*        1       
	 *        |  
	 *        |   
	 *   0 -- 2 -- 3       
	 */
	DistMat D1(4);
	D1(0,1) = 2; D1(0,2) = 1; D1(0,3) = 2;
	D1(1,2) = 1; D1(1,3) = 2;
	D1(2,3) = 1;
	int N = D1.size();
	DistMat W(D1, 2*N+2);
	Graph G1;
	std::vector<int> V({3});
	treerep_recurse(G1,W,V,0,1,2,N);
	G1.print();

	std::cout << "Test 2" << std::endl;
/*
*        1       
*      / | \
*     /  |  \
*   0 -- 2 -- 3       
*/
	DistMat D2(4);
	D2(0,1)=1; D2(0,2)=1; D2(0,3)=2;
	D2(1,2)=1; D2(1,3)=1;
	D2(2,3)=1;
	Graph G2;
	std::vector<int> V2({2});
	DistMat W2(D2, 2*N+2);
	treerep_recurse(G2,W2,V2,0,1,3,N);
	G2.print();
	return 0;
}
