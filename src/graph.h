#include <stdexcept>
#include <limits>
#include <algorithm>
#include <map>
#include <string>
#include <vector>
#include <random>
#include <fstream>
#include <iostream>

#ifndef GRAPH_H
#define GRAPH_H


/**
 * Symmetric matrix with 0's on the diagonal,
 * e.g. pairwise distances of N points. 
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
 * 		DistMat& operator *=(double d)
 * 			multiply all entries by a scalar
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
		int to_mtx(std::string file);
	private:
		unsigned _N;
		double _zero;
		std::vector<double> _data;
};

/**
 * An undirected graph stored in adjacency list. Vertices are labeled with ints.
 * 
 * Methods:
 * 		void add_edge(int u, int v)
 * 			Add edge (u,v). 
 * 		void remove_edge(int u, int v)
 *			Remove edge (u,v). 
 * 		void remove_vertex(int v)
 * 			Remove vertex v and all its edges. 
 * 		void retract(int u, int v)
 *			Retract edge (u,v) and label the new vertex u.
 *		std::vector<int> neighbors(int u)
 *			Get vector of vertices adjacent to u
 * 		DistMat metric()
 *			Calculate the shortest path distance between all vertices.
 *			Assumes vertices are 0...N and graph is connected
 *			@returns: symmetric matrix with distances
 * 		std::size_t size()
 * 			number of vertices in the graph	
 */
class Graph{
	public:
		void add_edge(int u, int v);
		void remove_edge(int u, int v);
		void remove_vertex(int v);
		void retract(int u, int v);
		std::vector<int> neighbors(int u);
		DistMat metric() const;
		std::size_t size() const;
		std::size_t num_edges() const;
		void print() const;
		int to_mtx(std::string file);
	private:
		void _rm(int u, int v);
		typedef std::map<int, std::vector<int> > vmap;
		typedef std::vector<int>::iterator vitr;
		typedef std::vector<int>::const_iterator const_vitr;
		vmap _adj;
};

/**
 * Generate a random tree
 * 		@param n: number of vertices in the tree
 */
Graph rand_tree(unsigned n);

 /* ======== MTX file utilities =========
 * For format spec see https://math.nist.gov/MatrixMarket/formats.html#MMformat 
 */
#define MTX_GRAPH_HDR "%MatrixMarket matrix coordinate pattern symmetric"
#define MTX_SYM_HDR "%MatrixMarket matrix coordinate real symmetric"

/**
 * int Graph::to_mtx(std::string file)
 * 		Write graph to mtx file as a `coordinate pattern symmetric' matrix
 * 		@param file: file name
 * 		@returns: 0 if successful, otherwise 1
 * 	int DistMat::to_mtx(std::string file)
 * 		Write matrix to mtx as `coordinate real symmetric' matrix
 * 		@returns: 0 if successful, otherwise 1 
 */

/** 
 * Load graph from mtx file
 * 		@param file: file with `coordinate pattern symmetric' matrix
 * 		@throws: runtime error if file is wrong format
 */
Graph graph_from_mtx(std::string file);

/**
 * Load DistMat from mtx file
 * 		@param file: file with `coordinate real symmetric' matrix 
 * 		@thows runtime error if file is wrong format 
 */
DistMat mat_from_mtx(std::string file);

#endif
