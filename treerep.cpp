#include "treerep.hpp"
#include <string>

double grmv_prod(int x, int y, int z, const DistMat& W){
	return 0.5*(W(x,z)+W(y,z)-W(x,y));
}

int treerep_recurse(Graph& G, DistMat& W, std::vector<int>& V, 
					int x, int y, int z, int N){
	static std::vector<int> stnr; // steiner nodes
	if (stnr.empty()){
		for (int i=2*N+2; i>=N; --i){
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
	if (std::abs(W(r,x)) < TOL && !retract){// retract (r, x)
		retract = true;
		W(r,x) = 0;
		G.remove_vertex(r);
		stnr.push_back(r);
		r = x;
		G.add_edge(x,y);
		G.add_edge(x,z);
	}

	W(r,y) = grmv_prod(y,x,z,W);
	if (std::abs(W(r,y)) < TOL &&  !retract){// retract (r, y)
		retract = true;
		W(r,x) = 0;
		W(r,y) = 0;
		G.remove_vertex(r);
		r = y ;
		stnr.push_back(r);

		G.add_edge(y,z);
		G.add_edge(y,x);
	}

	W(r,z) = grmv_prod(z,x,y,W);
	if (std::abs(W(r,y)) < TOL && !retract){
		//retract (r, z)
		retract = true;
		W(r,x) = 0;
		W(r,y) = 0;
		W(r,z) = 0;
		G.remove_vertex(r);
		stnr.push_back(r);
		r = z ;
		G.add_edge(y,z);
		G.add_edge(z,x);
	}

	// {zone1(r), zone1(z), zone2(z), zone1(y), zone2(y), zone1(x), zone2(x)}
	std::vector< std::vector<int> > zone(7);
	for (int i = 0; i< V.size(); ++i){
		int w = V[i];
		double a = grmv_prod(x,y,w,W); 
		double b = grmv_prod(y,z,w,W);
		double c = grmv_prod(z,x,w,W);
		int max = std::max({a,b,c});
		if ( std::abs(a-b) < TOL && std::abs(b-c) < TOL && std::abs(c-a) <TOL){
			if ( a < TOL &&  b < TOL && c<TOL && !retract){
				//retract (r,w)
				retract = true;
				G.remove_vertex(r);
				stnr.push_back(r);
				r = w;
				G.add_edge(w, x);
				G.add_edge(w, y);
				G.add_edge(w, z);
			}else{
				zone[0].push_back(w);	
				W(r,w) = (a+b+c)/3;
			}
		}else if (a == max){
			if (std::abs(W(z,w)-b) < TOL || std::abs(W(z,w)-c) < TOL){
				zone[1].push_back(w); // zone1(z)
			}else{
				zone[2].push_back(w); // zone2(z)
			}
			W(r,w) = a;
		}else if (b == max){
			if (std::abs(W(z,w)-a) < TOL || std::abs( W(z,w) - c) < TOL){
				zone[3].push_back(w); // zone1(y)
			}else{
				zone[4].push_back(w); // zone2(y)
			}
			W(r,w) = b;
		}else if (c == max){
			if (std::abs(W(z,w) - b) < TOL || std::abs(W(z,w) - a) < TOL){
				zone[5].push_back(w); // zone1(x)
			}else{
				zone[6].push_back(w); // zone2(x)
			}
			W(w,r) = c;
		}
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

void zone1_recurse(Graph& G, DistMat& W, std::vector<int>& V, int v, int N){
	int S = V.size(); 
	if (S == 1){
		int u = V.back();
		V.pop_back();
		G.add_edge(u,v);
	}else if (S >1){	
		int u = V.back();
		V.pop_back();
		int z = V.back();
		V.pop_back();
		treerep_recurse(G,W,V,v,u,z,N);
	}
}

void zone2_recurse(Graph& G, DistMat& W, std::vector<int>& V, int u, int v, int N){
	if (!V.empty()){
		int z = W.nearest(v, V);
		G.remove_edge(u,v);
		W(u,v) = 0;
		treerep_recurse(G,W,V,v,u,z,N);
	}
}

/* Remove vertex and all its edges.
 * Ignores vertices not in the graph*/
void Graph::remove_vertex(int v){
	vmap::iterator it = adj.find(v);
	if ( it != adj.end()){
		for (vmap::iterator gt=it; gt != adj.end(); ++gt){
			std::vector<int>::iterator j = std::find(gt->second.begin(), gt->second.end(), v);
			if (j != gt->second.end()){
				gt->second.erase(j);
			}
		}
		adj.erase(it);
	}
}

/* Add edge (u,v).
 * Note: only decreasing pairs u > v are stored in the adjacency list*/
void Graph::add_edge(int u, int v){ 
	if (u < v){
		add_edge(v,u);
	}else{
		std::vector<int>::iterator it = std::find(adj[u].begin(), adj[u].end(), v);
		if( it == adj[u].end()){
			adj[u].push_back(v);
		}
	}
}

/* Remove the edge (u,v) */
void Graph::remove_edge(int u, int v){
	if ( u < v){
		remove_edge(v,u);
	}else{
		std::vector<int>::iterator it = std::find(adj[u].begin(), adj[u].end(), v);
		if (it != adj[u].end()){
			adj[u].erase(it);
		}
	}
}

void Graph::print() const{
    for (vmap::const_iterator i=adj.begin(); i != adj.end(); ++i){
        std::cout << i->first << " :  ";
		for (std::vector<int>::const_iterator j = i->second.begin(); j != i->second.end(); ++j){
			std::cout << *j << " ";
		}
		std::cout << std::endl;
    }
    std::cout << "\n";
}

DistMat::DistMat(unsigned N): _N(N){
	_data = new double[(_N*(_N-1))/2];
	for (int i=0; i<_N; ++i){
		for (int j=i+1; j<_N; ++j){
			(*this)(i,j) = 0;
		}
	}
}

/* Inherit distances of M <= N points from another DistMat */
DistMat::DistMat(const DistMat& D, unsigned N): _N(N){
	int M = D.size();
	if (M > N){
		throw std::invalid_argument("Incompatible size");
	}
	_data = new double[(_N*(_N-1))/2];
	for (int i = 0; i<_N; ++i){
		for( int j = i+1; j<_N; ++j){
			if( i<M && j<M){
				(*this)(i,j) = D(i,j);
			}else{
				(*this)(i,j) = 0;
			}
		}
	}
}
	
 /* dist(x,y) = dist(y,x) and dist(x,x) = 0.
 * so we only need to store elements above the diagonal*/
double& DistMat::operator()(int i, int j){
	if(i >= _N || j >= _N || i < 0 || j < 0){
		throw std::invalid_argument("index out of bounds");
	}else if(i == j){
		throw std::invalid_argument("cannot assign to diagonal " + std::to_string(i));
	}else if (i > j){
		return _data[ _N*j + i - ((j+2)*(j+1))/2 ];
	}else{
		return _data[ _N*i + j - ((i+2)*(i+1))/2 ];
	}
}

double DistMat::operator()(int i, int j) const{
	if(i >= _N || j >= _N || i < 0 || j < 0){
		throw std::invalid_argument("index out of bounds");
	}else if(i == j){
		return 0;
	}else if (i >j){
		return _data[ _N*j + i - ((j+2)*(j+1))/2 ];
	}else{
		return _data[ _N*i + j - ((i+2)*(i+1))/2 ];
	}
}

/* find the point in pts nearest to i */
int DistMat::nearest(int i, const std::vector<int>& pts) const{
	if(pts.empty()){
		throw std::invalid_argument("set of points is empty");
	}else{
		double min = (*this)(i,pts.front());
		int n = 0;
		std::vector<int>::const_iterator it;
		for (it = pts.begin(); it != pts.end(); ++it){
			if( (*this)(i, *it) < min){
				min = (*this)(i, *it);
				n = *it;
			}
		}
		return n;
	}
}

int DistMat::size() const{
	return _N;
}

void DistMat::print() const{
	for (int i=0; i<_N; ++i){
		for (int j=0; j<_N; ++j){
			std::cout << (*this)(i,j) << " ";
		}
		std::cout << std::endl;
	}
}

DistMat::~DistMat(){
	delete[] _data;
}
