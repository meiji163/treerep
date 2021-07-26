#include <stdexcept>
#include <limits>
#include <algorithm>
#include <map>
#include <string>
#include <vector>
#include <set>
#include <iostream>

#ifndef GRAPH_H
#define GRAPH_H
/**
 * Represents a symmetric matrix with 0's on the diagonal,
 * e.g. a pairwise distances of N points. 
 * 
 * Constructors:
 * 		DistMat(unsigned N)
 * 			@param N: number of points ( > 0 )
 * 			@param val: value to fill matrix with (default=0)
 * 		DistMat(const DistMat& D, unsigned N) 
 * 			@param D: A size M <= N DistMat to copy values from
 * 		DistMat(const std::vector<double>& dist, unsigned N)
 * 			@param dist: A size N*(N-1)/2 vector of distances d(i,j), 0<=i<j<N
 * 						in lexographic order
 *
 * Methods:
 * 		double operator()(int i, int j)
 *			Access the element at (i,j). If i=j it is 0
 * 			@param i,j: ints >=0 and < N
 *		int nearest(int i, const std::vector<int>& pts) 
 *			Find the element of pts closest to point i. 
 *			@param pts: vector of non-negative ints < N
 *		int size()
 *			Return dimension of the matrix
 *		double max()
 *			Return max value in matrix	
 */
class DistMat{
	public:
		DistMat(unsigned N, double val=0);
		DistMat(const DistMat& D, unsigned N); 
		DistMat(const std::vector<double>& dist, unsigned N); 
		DistMat(double* dist, unsigned N);
		int nearest(int i, const std::vector<int>& pts) const; 
		double operator()(int i, int j) const;
		double& operator()(int i, int j);
		DistMat& operator*=(double d);
		double max() const;
		std::size_t size() const; 
		void print() const; 
	private:
		unsigned _N;
		double _zero;
		std::vector<double> _data;
};

/**
 * Undirected connected graph stored in adjacency list. 
 * Vertices are labeled with ints.
 */
class Graph{
	public:
		void add_edge(int u, int v);
		void remove_edge(int u, int v);
		void remove_vertex(int v);
		void retract(int u, int v);// Retract (u,v) and label the retracted vertex u.
		DistMat metric() const;
		std::size_t size() const;
		void print() const;
	private:
		void _rm(int u, int v);
		typedef std::map<int, std::vector<int> > vmap;
		typedef std::vector<int>::iterator vitr;
		typedef std::vector<int>::const_iterator const_vitr;
		vmap _adj;
};

/// load graph from mtx file
Graph load_graph(char* file);

/// TODO: save graph to mtx file
#endif
