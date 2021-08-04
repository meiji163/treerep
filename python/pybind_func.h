#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include "../src/graph.h"
#include "../src/treerep.h"

namespace py = pybind11;
typedef std::map<std::pair<int,int>, double> wmap;

py::array_t<double> py_metric(const Graph& G){
	DistMat D = G.metric();
	std::vector<double> metric = D.data();
	auto result = py::array_t<double>(metric.size());
	auto result_buf = result.request();
	int *result_ptr = (int *) result_buf.ptr;
	std::memcpy(result_ptr, metric.data(), metric.size()*sizeof(double));
	return result;
}

std::pair<Graph, wmap> py_graph_treerep(const Graph& G, double tol=0.1){
	DistMat D = G.metric();
	int N = D.size();
	std::pair<Graph,DistMat> pr = treerep(G.metric(),tol);
	Graph T = pr.first;
	DistMat W = pr.second;
	wmap weight;
	for(int i=0; i<W.size(); ++i){
		for( int j=i+1; j<W.size(); ++j){
			if(T.is_adj(i,j)){
				weight[std::make_pair(i,j)]=W(i,j);
			}
		}
	}
	return std::make_pair(T,weight);
}

// accepts 1D numpy array of size N*(N-1)/2 as the N-point metric 
std::pair<Graph, wmap> py_treerep(py::array_t<double, py::array::c_style | py::array::forcecast> metric, int N, double tol=0.1){
	if( N<2 || metric.size() != (N*(N-1))/2){
		throw std::invalid_argument("array must be size N*(N-1)/2");
	}
	DistMat M(metric.data(), N);
	std::pair<Graph,DistMat> pr = treerep(M,tol);
	Graph T = pr.first;
	DistMat W = pr.second;
	wmap weight;
	for(int i=0; i<W.size(); ++i){
		for( int j=i+1; j<W.size(); ++j){
			if(W(i,j)<0){
				W(i,j)=0;
			}
			if(T.is_adj(i,j)){
				weight[std::make_pair(i,j)]=W(i,j);
			}
		}
	}

	return std::make_pair(T,weight);
}

