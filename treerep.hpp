#include <stdexcept>
#include <algorithm>
#include <cmath>
#include <map>
#include <vector>
#include <iostream>
#define TOL 0.07

/* represents a undirected graph */
struct Graph{
	void add_edge(int u, int v);
	void remove_edge(int u, int v);
	void remove_vertex(int v);
	void print() const;
	typedef std::map<int, std::vector<int> > vmap;
	vmap adj;
};

/* Holds pairwise distances of N points */
class DistMat{
	public:
		DistMat(unsigned N);
		DistMat(const DistMat& D, unsigned N);
		double operator()(int i, int j) const;
		double& operator()(int i, int j);
		int nearest(int i, const std::vector<int>& pts) const;
		int size() const;
		void print() const;
		~DistMat();
	private:
		unsigned _N;
		double* _data;
};

int treerep_recurse(Graph& G, DistMat& W, std::vector<int>& V, int x, int y, int z, int N);
void zone1_recurse(Graph& G, DistMat& W, std::vector<int>& V, int v, int N);
void zone2_recurse(Graph& G, DistMat& W, std::vector<int>& V, int u, int v, int N);
double grmv_prod(int x, int y, int z, const DistMat& W);
