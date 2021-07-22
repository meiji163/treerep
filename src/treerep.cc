#include <algorithm>
#include <map>
#include <cmath>
#include <vector>
#include <set>
#include <random>
#include <iostream>
#include "graph.h"

#define TR_NUM_THREADS 4
std::random_device rd;
std::default_random_engine rng(rd());

int treerep_recurse(Graph& G, DistMat& W, std::vector<int>& V, 
					int x, int y, int z, int N, double tol=0.1);
void zone1_recurse(Graph& G, DistMat& W, std::vector<int>& V, int v, 
					int N, double tol=0.1);
void zone2_recurse(Graph& G, DistMat& W, std::vector<int>& V, int u, int v, 
					int N, double tol=0.1);
void sort(Graph& G, DistMat& W, std::vector<int>& V, std::vector< std::vector<int> >& zone, int x, int y, int z, 
		  int r, int N, bool retract, double tol, std::vector<int>& stnr);

void print(const std::vector<int>& vec){
	for (std::vector<int>::const_iterator i = vec.begin(); i!= vec.end(); ++i){
		std::cout << *i << " ";
	}
	std::cout << std::endl;
}

double grmv_prod(int x, int y, int z, const DistMat& W){
	return 0.5*(W(x,z)+W(y,z)-W(x,y));
}

int treerep_recurse(Graph& G, DistMat& W, std::vector<int>& V, 
					int x, int y, int z, int N, double tol){
	static std::vector<int> stnr; // steiner nodes
	if (stnr.empty()){
		for (int i=W.size(); i>=N; --i){
			stnr.push_back(i);
		}
	}

	// form the universal tree for x,y,z
	int r = stnr.back();
	stnr.pop_back();
	G.add_edge(r,x);
	G.add_edge(r,y);
	G.add_edge(r,z);

	bool retract = false;
	W(r,x) = grmv_prod(x,y,z,W); 
	if (std::abs(W(r,x)) < tol && !retract){// retract (r, x)
		W(r,x) = 0;
		G.retract(x,r);
		stnr.push_back(r);
		r = x;
		retract = true;
	}
	W(r,y) = grmv_prod(y,x,z,W);
	if (std::abs(W(r,y)) < tol &&  !retract){// retract (r, y)
		W(r,x) = 0;
		W(r,y) = 0;
		G.retract(y,r);
		stnr.push_back(r);
		r = y;
		retract = true;
	}
	W(r,z) = grmv_prod(z,x,y,W);
	if (std::abs(W(r,z)) < tol && !retract){ //retract (r, z)
		W(r,x) = 0;
		W(r,y) = 0;
		W(r,z) = 0;
		G.retract(z,r);
		stnr.push_back(r);
		r = z ;
		retract = true;
	}

	// {zone1(r), zone1(z), zone2(z), zone1(y), zone2(y), zone1(x), zone2(x)}
	std::vector< std::vector<int> > zone(7);
	sort(G,W,V,zone,x,y,z,r,
			N,retract,tol,stnr);

	for (int i=0; i<7; ++i){
		std::cout << "zone" << i << ": " << std::endl;
		print(zone[i]);
	}
	zone1_recurse(G,W,zone[0],r,N);
	zone1_recurse(G,W,zone[1],z,N);
	zone1_recurse(G,W,zone[3],y,N);
	zone1_recurse(G,W,zone[5],x,N);
	zone2_recurse(G,W,zone[2],z,r,N);
	zone2_recurse(G,W,zone[4],y,r,N);
	zone2_recurse(G,W,zone[6],x,r,N);
	return 0;
}

void sort(Graph& G, DistMat& W, std::vector<int>& V, std::vector< std::vector<int> >& zone, 
		int x, int y, int z, int r, int N, bool retract, double tol, std::vector<int>& stnr){
	for (int i = 0; i< V.size(); ++i){
		//std::cout << "(r,x,y,z) = " << r << " " << x << " " << y << " " << z << std::endl;
		int w = V[i];
		//std::cout << "w = " << w << std::endl;
		double a = grmv_prod(x,y,w,W); 
		double b = grmv_prod(y,z,w,W);
		double c = grmv_prod(z,x,w,W);
		double max = std::max({a,b,c});
		if ( std::abs(a-b) < tol && std::abs(b-c) < tol && std::abs(c-a) <tol){
			if ( a < tol &&  b < tol && c<tol && !retract){ //retract (r,w)
				retract = true;
				for( int i = N; i<W.size(); ++i){
					W(w,i) = W(r,i);
				}
				W(r,x) = 0;
				W(r,y) = 0;
				W(r,z) = 0;
				G.retract(w,r);
				stnr.push_back(r);
				r = w;
			}else{
				zone[0].push_back(w); //zone1(r)
				W(r,w) = (a+b+c)/3;
			}
		}else if (a == max){
			if (std::abs(W(z,w)-b) < tol || std::abs(W(z,w)-c) < tol){
				zone[1].push_back(w); // zone1(z)
			}else{
				zone[2].push_back(w); // zone2(z)
			}
			W(r,w) = a; 
		}else if (b == max){
			if (std::abs(W(z,w)-a) < tol || std::abs( W(z,w) - c) < tol){
				zone[3].push_back(w); // zone1(y)
			}else{
				zone[4].push_back(w); // zone2(y)
			}
			W(r,w) = b;
		}else if (c == max){
			if (std::abs(W(z,w) - b) < tol || std::abs(W(z,w) - a) < tol){
				zone[5].push_back(w); // zone1(x)
			}else{
				zone[6].push_back(w); // zone2(x)
			}
			W(w,r) = c;
		}
	}
}

void zone1_recurse(Graph& G, DistMat& W, std::vector<int>& V, int v, 
					int N, double tol){
	int S = V.size(); 
	if (S == 1){
		int u = V.back();
		V.pop_back();
		G.add_edge(u,v);
	}else if (S >1){	
		std::shuffle(V.begin(), V.end(), rng);
		int u = V.back();
		V.pop_back();
		int z = V.back();
		V.pop_back();
		treerep_recurse(G,W,V,v,u,z,N);
	}
}

void zone2_recurse(Graph& G, DistMat& W, std::vector<int>& V, int u, int v, 
					int N, double tol){
	if (!V.empty()){
		int z = W.nearest(v, V);
		G.remove_edge(u,v);
		W(u,v) = 0;
		treerep_recurse(G,W,V,v,u,z,N);
	}
}

int test1(){
	std::cout << "Test 1:" << std::endl;
/*        1       
 *        |  
 *        |   
 *   0 -- 2 -- 3       
 */
	int N = 4; 
	DistMat D({2,1,2,
			   1,2,
			   1}, N);
	DistMat W(D, 2*N+2);
	Graph G;
	std::vector<int> V({3});
	treerep_recurse(G,W,V,0,1,2,N);
	G.print();
	return 0;
}

int test2(){
	std::cout << "Test 2:" << std::endl;
/*
*        1       
*      / | \
*     /  |  \
*   0 -- 2 -- 3       
*/
	int N =4;
	DistMat D({1,1,2,
			   1,1,
			   1}, N);
	Graph G;
	std::vector<int> V({2});
	DistMat W(D, 2*N+2);
	treerep_recurse(G,W,V,0,1,3,N);
	G.print();
	return 0;
}

int test3(){
	std::cout << "Test 3:" << std::endl;
	int N = 8;
	DistMat D( {32,48,51,50,48,98,148,
				26,34,29,33,84,136,
				42,44,44,92,152,
				44,38,86,142,
				42,89,142,
				90,142,
				148},
				N);
	DistMat W(D, 3*N);
	double m = W.max();
	W *= 1/m;
	Graph G;
	std::vector<int> V({1,3,5,6,7});
	treerep_recurse(G,W,V,0,2,4,N);
	G.print();
	return 0;
}

int main(){
	return test3();
}

