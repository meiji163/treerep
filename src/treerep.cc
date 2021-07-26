#include "treerep.h"
unsigned TR_N;
double TR_TOL;
std::default_random_engine TR_RNG(1);

void print(const std::vector<int>& vec){
	for (std::vector<int>::const_iterator i = vec.begin(); i!= vec.end(); ++i){
		std::cout << *i << " ";
	}
	std::cout << std::endl;
}

std::pair<Graph,DistMat> treerep(const DistMat& D, double tol){
	TR_TOL = tol;
	TR_N = D.size();
	//normalize metric 
	double max = D.max();
	DistMat W(D, 2*TR_N);
	W *= 1/max;
	std::vector<int> V(TR_N);
	for (int i=0; i<TR_N; ++i){
		V[i] = i;
	}

	//choose random starting vertices
	std::shuffle(V.begin(), V.end(), TR_RNG);
	int x = V.back();
	V.pop_back();
	int y = V.back();
	V.pop_back();
	int z = V.back();
	V.pop_back();

	//initialize steiner nodes
	std::vector<int> stn;
	stn.reserve(TR_N);
	for (int i=2*TR_N; i>=TR_N; --i){
		stn.push_back(i);
	}
	Graph G;
	int ex = _treerep_recurse(G,W,V,stn,x,y,z);
	if(ex){
		std::cout << "Error code " << ex << std::endl;
	}else{
		W *= max;
	}
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
	if (std::abs(W(r,x)) < TR_TOL && !rtr){// retract (r, x)
		W(r,x) = 0;
		G.retract(x,r);
		stn.push_back(r);
		r = x;
		rtr = true;
	}
	W(r,y) = grmv_prod(y,x,z,W);
	if (std::abs(W(r,y)) < TR_TOL && !rtr){// retract (r, y)
		W(r,x) = 0;
		W(r,y) = 0;
		G.retract(y,r);
		stn.push_back(r);
		r = y;
		rtr = true;
	}
	W(r,z) = grmv_prod(z,x,y,W);
	if (std::abs(W(r,z)) < TR_TOL && !rtr){ //retract (r, z)
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

	//for (int i=0; i<7; ++i){
		//std::cout << "zone" << i << ": " << std::endl;
		//print(zone[i]);
	//}
	_zone1_recurse(G,W,zone[0],stn,r);
	_zone1_recurse(G,W,zone[1],stn,z);
	_zone1_recurse(G,W,zone[3],stn,y);
	_zone1_recurse(G,W,zone[5],stn,x);
	_zone2_recurse(G,W,zone[2],stn,z,r);
	_zone2_recurse(G,W,zone[4],stn,y,r);
	_zone2_recurse(G,W,zone[6],stn,x,r);
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
		if ( std::abs(a-b)<TR_TOL && std::abs(b-c)<TR_TOL && std::abs(c-a)<TR_TOL){
			if (a<TR_TOL &&  b<TR_TOL && c<TR_TOL && !rtr){ //retract (r,w)

				std::cout << "retract " << r << " " << w << std::endl;

				rtr = true;
				for( int i = TR_N; i<W.size(); ++i){
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
			if (std::abs(W(z,w)-b) < TR_TOL || std::abs(W(z,w)-c) < TR_TOL){
				zone[1].push_back(w); // zone1(z)
			}else{
				zone[2].push_back(w); // zone2(z)
			}
			W(r,w) = a; 
		}else if (b == max){
			if (std::abs(W(z,w)-a) < TR_TOL || std::abs( W(z,w) - c) < TR_TOL){
				zone[3].push_back(w); // zone1(y)
			}else{
				zone[4].push_back(w); // zone2(y)
			}
			W(r,w) = b;
		}else if (c == max){
			if (std::abs(W(z,w) - b) < TR_TOL || std::abs(W(z,w) - a) < TR_TOL){
				zone[5].push_back(w); // zone1(x)
			}else{
				zone[6].push_back(w); // zone2(x)
			}
			W(w,r) = c;
		}
	}
	return zone;
}

void _zone1_recurse(Graph& G, DistMat& W, std::vector<int>& V, std::vector<int>& stn, int v){
	int S = V.size(); 
	if (S == 1){
		int u = V.back();
		V.pop_back();
		G.add_edge(u,v);
	}else if (S>1){	
		std::shuffle(V.begin(), V.end(), TR_RNG);
		int u = V.back();
		V.pop_back();
		int z = V.back();
		V.pop_back();
		_treerep_recurse(G,W,V,stn,v,u,z);
	}
}

void _zone2_recurse(Graph& G, DistMat& W, std::vector<int>& V, std::vector<int>& stn, 
					int u, int v){
	if (!V.empty()){
		int z = W.nearest(v, V);
		G.remove_edge(u,v);
		W(u,v) = 0;
		_treerep_recurse(G,W,V,stn,v,u,z);
	}
}

inline double grmv_prod(int x, int y, int z, const DistMat& W){
	return 0.5*(W(x,z)+W(y,z)-W(x,y));
}
