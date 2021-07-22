#include "graph.h"

void Graph::remove_vertex(int v){
	vmap::iterator i = adj.find(v);
	if (i != adj.end()){
		std::set<int>::iterator j;
		for (j = adj[v].begin(); j!= adj[v].end(); ++j){
			adj[*j].erase(v);
		}
		adj.erase(i);
	}
}

void Graph::add_edge(int u, int v){ 
	std::set<int>::iterator it = std::find(adj[u].begin(), adj[u].end(), v);
	if( it == adj[u].end()){
		adj[u].insert(v);
		adj[v].insert(u);
	}
}

void Graph::remove_edge(int u, int v){
	std::set<int>::iterator it = std::find(adj[u].begin(), adj[u].end(), v);
	if (it != adj[u].end()){
		adj[u].erase(it);
		it = std::find(adj[v].begin(), adj[v].end(), u);
		if (it != adj[v].end()){
			adj[v].erase(it);
		}
	}
}

void Graph::retract(int u, int v){
	std::set<int>::iterator it;
	for ( it = adj[v].begin(); it!=adj[v].end(); ++it){
		adj[u].insert(*it);
		adj[*it].insert(u);
		adj[*it].erase(v);
	}
	adj[u].erase(v);
	adj.erase(v);
}

void Graph::print() const{
    for (vmap::const_iterator i=adj.begin(); i != adj.end(); ++i){
        std::cout << i->first << " :  ";
		std::set<int>::const_iterator j;
		for (j = i->second.begin(); j != i->second.end(); ++j){
			std::cout << *j << " ";
		}
		std::cout << std::endl;
    }
    std::cout << "\n";
}

DistMat::DistMat(unsigned N): _N(N), _zero(0){
	_data = new double[(_N*(_N-1))/2];
	for (int i=0; i<_N; ++i){
		for (int j=i+1; j<_N; ++j){
			(*this)(i,j) = 0;
		}
	}
}

DistMat::DistMat(const DistMat& D, unsigned N):_N(N), _zero(0.){
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

DistMat::DistMat(const std::vector<double>& dist, unsigned N): _N(N), _zero(0.){
	int S = dist.size();
	if (S != (N*(N-1))/2){
		throw std::invalid_argument("Incompatible sizes "+std::to_string(S) 
									+" and "+std::to_string(N));
	}
	_data = new double[(_N*(_N-1))/2];
	for (int i=0; i<S; ++i){
		_data[i] = dist[i];
	}
}

double& DistMat::operator()(int i, int j){
	if(i >= _N || j >= _N || i < 0 || j < 0){
		throw std::invalid_argument("index out of bounds");
	}else if(i == j){
		return _zero;
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

DistMat& DistMat::operator*=(double d){
	for(int i=0; i< (_N*(_N-1))/2; ++i){
		_data[i] *= d;
	}
	return *this;
}

double DistMat::max(){
	return *std::max_element(_data, _data+ (_N*(_N-1))/2);
}

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
