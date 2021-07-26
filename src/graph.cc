#include "graph.h"

void Graph::remove_vertex(int v){
	vmap::iterator i = _adj.find(v);
	if (i != _adj.end()){
		vitr j;
		for (j = _adj[v].begin(); j!= _adj[v].end(); ++j){
			_rm(*j,v);
		}
		_adj.erase(i);
	}
}

void Graph::add_edge(int u, int v){ 
	vitr it = std::find(_adj[u].begin(), _adj[u].end(), v);
	if( it == _adj[u].end()){
		_adj[u].push_back(v);
		_adj[v].push_back(u);
	}
}

void Graph::remove_edge(int u, int v){
	_rm(u,v);
	_rm(v,u);
}

void Graph::retract(int u, int v){
	for (vitr it = _adj[v].begin(); it!=_adj[v].end(); ++it){
		if( u!=*it){
			_adj[u].push_back(*it);
			_rm(*it,v);
			_adj[*it].push_back(u);
		}
	}
	_rm(u,v);
	_adj.erase(v);
}

inline void Graph::_rm(int u, int v){
	vitr it = std::find(_adj[u].begin(), _adj[u].end(), v);
	if (it != _adj[u].end()){
		_adj[u].erase(it);
	}
}

void Graph::print() const{
	std::cout << "======= GRAPH =======" << std::endl;
    for (vmap::const_iterator i=_adj.begin(); i != _adj.end(); ++i){
        std::cout << i->first << " :  ";
		const_vitr j;
		for (j = i->second.begin(); j != i->second.end(); ++j){
			std::cout << *j << " ";
		}
		std::cout << std::endl;
    }
    std::cout << "\n";
}

inline std::vector<int> Graph::neighbors(int u){
	return _adj[u];
}

DistMat Graph::metric() const{
	double infty = std::numeric_limits<double>::infinity();
	DistMat W(_adj.size(), infty); 
	vmap::const_iterator itr,jtr,ktr;
	for (itr=_adj.begin(); itr!=_adj.end();++itr){
		for(const_vitr vit=itr->second.begin();vit!=itr->second.end(); ++vit){ 
			W(itr->first,*vit) = 1;
		}
	}

	// Floyd-Warshall
	int i,j,k;
	for (k=0,ktr=_adj.begin(); ktr!=_adj.end(); ++ktr,++k){
		for (i=0,itr=_adj.begin(); itr!=_adj.end(); ++itr,++i){
			for (j=i+1,jtr=_adj.begin(); j<_adj.size();++jtr,++j){
				if( i!=j && j!=k && k!=i
					&& (W(i,j) > W(i,k) + W(k,j)) ){
					W(i,j) = W(i,k) + W(k,j);
				}
			}
		}
	}
	return W;
}

std::size_t Graph::size() const{
	return _adj.size();
}

std::size_t Graph::num_edges() const{
	unsigned len=0;
	for (vmap::const_iterator it=_adj.begin();it!=_adj.end(); ++it){
		len += (it->second).size();
	}
	return len/2;
}

Graph graph_from_mtx(std::string file){
	Graph G;
	std::ifstream f(file, std::ios::in);
	std::string s;
	if (!f.is_open()){
		throw std::runtime_error(std::string() + "Can't open file: " + file);
	}
	std::getline(f,s);
	if (s.find(MTX_GRAPH_HDR) == std::string::npos){
		throw std::runtime_error(std::string() + "Invalid mtx file: "+ file);
	}
	//skip comments
	char c = f.peek();
	while(c=='%'){
		std::getline(f,s);
		c = f.peek();
	}
	int rows, cols, lns;
	f >> rows >> cols >> lns;
	if (rows != cols){
		throw std::runtime_error(std::string()+"Rows not equal to columns: " + file);
	}
	int u,v;
	for (int i=0; i<lns; ++i){
		f >> u >> v;
		G.add_edge(u-1,v-1); // 1-indexed labels 
	}
	f.close();
	return G;
}

DistMat mat_from_mtx(std::string file){
	std::ifstream f(file, std::ios::in);
	std::string s;
	if (!f.is_open()){
		throw std::runtime_error(std::string() + "Can't open file: " + file);
	}
	std::getline(f,s);
	if (s.find(MTX_SYM_HDR)==std::string::npos){
		throw std::runtime_error(std::string() + "Invalid mtx file: "+ file);
	}
	//skip comments
	char c = f.peek();
	while(c=='%'){
		std::getline(f,s);
		c = f.peek();
	}
	//rows, cols, entries
	int rows, cols, lns;
	f >> rows >> cols >> lns;
	if (rows != cols){
		throw std::runtime_error(std::string()+"Rows not equal to columns: " + file);
	}

	DistMat M(rows);
	int u,v;
	double val;
	for(int i=0; i<lns; ++i){
		f >> u >> v >> val;
		M(u-1,v-1) = val;
	}
	f.close();
	return M;
}


int Graph::to_mtx(std::string file){
	std::ofstream f(file,std::ios::out);
	if(!f.is_open()){
		return 1;
	}
	f << MTX_GRAPH_HDR << std::endl;
	int S = this->size();
	f << S << " " << S << " " << this->num_edges() << std::endl;
	for (int u=0; u<_adj.size(); ++u){
		std::vector<int> nbr = neighbors(u);
		for(vitr vt=nbr.begin(); vt!=nbr.end(); ++vt){
			if(*vt > u){
				f << u+1 << " " << (*vt)+1 << std::endl;
			}
		}
	}
	f.close();
	return 0;
}

int DistMat::to_mtx(std::string file){ 
	std::ofstream f(file, std::ios::out);
	if (!f.is_open()){
		return 1;
	}
	f << MTX_SYM_HDR << std::endl;
	f << _N << " " << _N << " " << (_N*(_N-1))/2 << std::endl;
	for (int i=0; i<_N; ++i){
		for (int j=i+1; j<_N; ++j){
			f << i+1 << " " << j+1 << " " << (*this)(i,j) << std::endl;
		}
	}
	return 0;
}

DistMat::DistMat(unsigned N, double val): _N(N), _zero(0){
	_data.resize((N*(N-1))/2);
	for (int i=0; i<_N; ++i){
		for (int j=i+1; j<_N; ++j){
			(*this)(i,j) = val;
		}
	}
}

DistMat::DistMat(const DistMat& D, unsigned N):_N(N), _zero(0.){
	int M = D.size();
	if (M > N){
		throw std::invalid_argument("Incompatible size");
	}
	_data.resize((N*(N-1))/2);
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
	_data.resize((N*(N-1))/2);
	for (int i=0; i<S; ++i){
		_data[i] = dist[i];
	}
}

DistMat::DistMat(double* dist, unsigned N): _N(N), _zero(0.){
	_data.resize((N*(N-1))/2);
	for (int i=0; i< (N*(N-1))/2; ++i){
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

double DistMat::max() const{
	return *std::max_element(_data.begin(), _data.end());
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

inline std::size_t DistMat::size() const{
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
