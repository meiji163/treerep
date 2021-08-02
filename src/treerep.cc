#include "treerep.h"
unsigned TREP_N;
double TREP_TOL;

std::default_random_engine& _trep_rng(int seed){
	static std::default_random_engine rng(seed);
	return rng;
}

std::pair<Graph,DistMat> treerep(const DistMat& D, double tol){
	TREP_N = D.size();
	double max = D.max();
	TREP_TOL = tol*max;
	DistMat W(D, 2*TREP_N);
	std::vector<int> V(TREP_N);
	for (int i=0; i<TREP_N; ++i){
		V[i] = i;
	}

	//choose random starting vertices
	std::shuffle(V.begin(), V.end(), _trep_rng());
	int x = V.back();
	V.pop_back();
	int y = V.back();
	V.pop_back();
	int z = V.back();
	V.pop_back();

	//initialize steiner nodes
	std::vector<int> stn;
	stn.reserve(TREP_N);
	for (int i=2*TREP_N; i>=TREP_N; --i){
		stn.push_back(i);
	}
	Graph G;
	_treerep_recurse(G,W,V,stn,x,y,z);
	std::pair<Graph,DistMat> pr = std::make_pair(G,W);
	return pr;
}

int _treerep_recurse(Graph& G, DistMat& W, std::vector<int>& V, std::vector<int>& stn,
					int x, int y, int z){
	if(stn.empty()){
		return 1;
	}

	// form the universal tree for x,y,z
	int r = stn.back();
	stn.pop_back();
	G.add_edge(r,x);
	G.add_edge(r,y);
	G.add_edge(r,z);

	bool rtr = false;
	W(r,x) = grmv_prod(x,y,z,W); 
	if (std::abs(W(r,x)) < TREP_TOL && !rtr){// retract (r, x)
		W(r,x) = 0;
		G.retract(x,r);
		stn.push_back(r);
		r = x;
		rtr = true;
	}
	W(r,y) = grmv_prod(y,x,z,W);
	if (std::abs(W(r,y)) < TREP_TOL && !rtr){// retract (r, y)
		W(r,x) = 0;
		W(r,y) = 0;
		G.retract(y,r);
		stn.push_back(r);
		r = y;
		rtr = true;
	}
	W(r,z) = grmv_prod(z,x,y,W);
	if (std::abs(W(r,z)) < TREP_TOL && !rtr){ //retract (r, z)
		W(r,x) = 0;
		W(r,y) = 0;
		W(r,z) = 0;
		G.retract(z,r);
		stn.push_back(r);
		r = z ;
		rtr = true;
	}
	//sort rest of vertices into 7 zones
	vecvec zone = _sort(G,W,V,stn,x,y,z,r,rtr);

	_zone1(G,W,zone[0],stn,r);
	_zone1(G,W,zone[1],stn,z);
	_zone1(G,W,zone[3],stn,y);
	_zone1(G,W,zone[5],stn,x);
	_zone2(G,W,zone[2],stn,z,r);
	_zone2(G,W,zone[4],stn,y,r);
	_zone2(G,W,zone[6],stn,x,r);
	return 0;
}

vecvec _sort(Graph& G, DistMat& W, std::vector<int>& V, std::vector<int>& stn,
			int x, int y, int z, int r, bool rtr){ 
	vecvec zone(7);
	for (int i = 0; i< V.size(); ++i){
		int w = V[i];
		double a = grmv_prod(x,y,w,W); 
		double b = grmv_prod(y,z,w,W);
		double c = grmv_prod(z,x,w,W);
		double max = std::max({a,b,c});
		if ( std::abs(a-b)<TREP_TOL && std::abs(b-c)<TREP_TOL && std::abs(c-a)<TREP_TOL){
			if (a<TREP_TOL &&  b<TREP_TOL && c<TREP_TOL && !rtr){ //retract (r,w)
				rtr = true;
				for( int i = TREP_N; i<W.size(); ++i){
					W(w,i) = W(r,i);
				}
				W(r,x) = 0;
				W(r,y) = 0;
				W(r,z) = 0;
				G.retract(w,r);
				stn.push_back(r);
				r = w;
			}else{
				zone[0].push_back(w); //zone1(r)
				W(r,w) = (a+b+c)/3;
			}
		}else if (a == max){
			if (std::abs(W(z,w)-b)<TREP_TOL || std::abs(W(z,w)-c)<TREP_TOL){
				zone[1].push_back(w); // zone1(z)
			}else{
				zone[2].push_back(w); // zone2(z)
			}
			W(r,w) = a; 
		}else if (b == max){
			if (std::abs(W(z,w)-a)<TREP_TOL || std::abs(W(z,w)-c)<TREP_TOL){
				zone[3].push_back(w); // zone1(y)
			}else{
				zone[4].push_back(w); // zone2(y)
			}
			W(r,w) = b;
		}else if (c == max){
			if (std::abs(W(z,w)-b)<TREP_TOL || std::abs(W(z,w)-a)<TREP_TOL){
				zone[5].push_back(w); // zone1(x)
			}else{
				zone[6].push_back(w); // zone2(x)
			}
			W(w,r) = c;
		}
	}
	return zone;
}

void _zone1(Graph& G, DistMat& W, std::vector<int>& V, std::vector<int>& stn, int v){
	int S = V.size(); 
	if (S == 1){
		int u = V.back();
		V.pop_back();
		G.add_edge(u,v);
	}else if (S>1){	
		std::shuffle(V.begin(), V.end(), _trep_rng());
		int u = V.back();
		V.pop_back();
		int z = V.back();
		V.pop_back();
		_treerep_recurse(G,W,V,stn,v,u,z);
	}
}

void _zone2(Graph& G, DistMat& W, std::vector<int>& V, std::vector<int>& stn, 
					int u, int v){
	if (!V.empty()){
		int z = W.nearest(v, V);
		G.remove_edge(u,v);
		_treerep_recurse(G,W,V,stn,v,u,z);
	}
}

inline double grmv_prod(int x, int y, int z, const DistMat& W){
	return 0.5*(W(x,z)+W(y,z)-W(x,y));
}

Graph rand_tree(unsigned n, int seed){
	std::uniform_int_distribution<> dist(0,n);
	Graph G;
	int v;
	for (int u=1; u<n; ++u){
		v = dist(_trep_rng()) %u;
		G.add_edge(u,v);
	}
	return G;
}

