#include <stdexcept>
#include <algorithm>
#include <map>
#include <string>
#include <vector>
#include <set>
#include <iostream>

#pragma once
/**
 * Undirected graph stored in adjacency list. Vertices are labeled with ints.
 * Insert/remove functions do nothing if vertex/edge is already in graph (resp. not in graph)
 */
struct Graph{
	void add_edge(int u, int v);// Insert the edge (u,v)
	void remove_edge(int u, int v);// Remove the edge (u,v) 
	void remove_vertex(int v);// Remove vertex v 
	void retract(int u, int v);// Retract (u,v) and label the retracted vertex u.
	void print() const;
	typedef std::map<int, std::set<int> > vmap;
	vmap adj;
};

/**
 * Represents matrix of pairwise distances for N points. 
 * In other words, a symmetric matrix with 0's on the main diagonal.
 * 
 * Constructors:
 * 		DistMat(unsigned N)
 * 			@param N: number of points ( > 0 )
 * 		DistMat(const DistMat& D, unsigned N) 
 * 			@param D: A size M <= N DistMat to copy values from
 * 		DistMat(const std::vector<double>& dist, unsigned N)
 * 			@param dist: A size N*(N-1)/2 vector of distances d(i,j), 0<=i<j<N
 * 						in lexographic order
 *
 * Methods:
 * 		double operator()(int i, int j)
 * 			@param i,j: integers >=0 and < N
 *			access the element at (i,j). If i=j it is 0
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
		DistMat(unsigned N);
		DistMat(const DistMat& D, unsigned N); 
		DistMat(const std::vector<double>& dist, unsigned N); 
		int nearest(int i, const std::vector<int>& pts) const; 
		double operator()(int i, int j) const;
		double& operator()(int i, int j);
		DistMat& operator*=(double d);
		double max();
		int size() const; 
		void print() const; 
		~DistMat();
	private:
		unsigned _N;
		double* _data;
		double _zero;
};
