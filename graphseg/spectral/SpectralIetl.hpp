/*
 * SpectralIetl.hpp
 *
 *  Created on: Jan 19, 2012
 *      Author: david
 */

#ifndef DASP_SPECTRAL_SPECTRALIETL_HPP_
#define DASP_SPECTRAL_SPECTRALIETL_HPP_

#include "../Common.hpp"
#include "../as_range.hpp"
#include <boost/graph/adjacency_list.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/io.hpp>
// #include <ietl/interface/ublas.h>
// #include <ietl/vectorspace.h>
// #include <ietl/lanczos.h>
#include <boost/random.hpp>
#include <boost/limits.hpp>
#include <limits>
#include <iostream>
#include <vector>
#include <cmath>

namespace graphseg { namespace detail {

template<typename Graph>
std::vector<EigenComponent> SpectralIetl(const Graph& graph, unsigned int num_ev)
{
#if 0
#ifdef SPECTRAL_VERBOSE
	std::cout << "Sparse Solver: started" << std::endl;
#endif

	// We want to solve the EV problem: (D - W) x = \lamda D x.
	// Each edge of the graph defines two entries into the symmetric matrix W.
	// The diagonal matrix D is defined via d_i = sum_j{w_ij}.

	// As D is a diagonal matrix the the general problem can be easily transformed
	// into a normal eigenvalue problem by decomposing D = L L^t, which yields L = sqrt(D).
	// Thus the EV problem is: L^{-1} (D - W) L^{-T} y = \lambda y.
	// Eigenvectors can be transformed using x = L^{-T} y.

#ifdef SPECTRAL_VERBOSE
	std::cout << "Sparse Solver: preparing problem" << std::endl;
#endif

	// The dimension of the problem
	const int n = boost::num_vertices(graph);

	struct Entry {
		int i, j;
		Real value;
	};

	// Each edge defines two entries (one in the upper and one in the lower).
	// In addition all diagonal entries are non-zero.
	// Thus the number of non-zero entries in the lower triangle is equal to
	// the number of edges plus the number of nodes.
	// This is not entirely true as some connections are possibly rejected.
	// Additionally some connections may be added to assure global connectivity.
	const int nnz_guess = boost::num_edges(graph) + n;

	// collect all non-zero elements
	std::vector<Entry> entries;
	entries.reserve(nnz_guess);

	// also collect diagonal entries
	std::vector<Real> diag(n);

	// no collect entries
	for(auto eid : as_range(boost::edges(graph))) {
		int ea = static_cast<int>(boost::source(eid, graph));
		int eb = static_cast<int>(boost::target(eid, graph));
		float ew = boost::get(boost::edge_weight, graph, eid);
		// assure correct edge weight
		if(std::isnan(ew)) {
			std::cerr << "Weight for edge (" << ea << "," << eb << ") is nan!" << std::endl;
			continue;
		}
		if(ew < 0) {
			std::cerr << "Weight for edge (" << ea << "," << eb << ") is negative!" << std::endl;
			continue;
		}
		// assure that no vertices is connected to self
		if(ea == eb) {
			std::cerr << "Vertex " << ea << " is connected to self!" << std::endl;
			continue;
		}
		// In the lower triangle the row index i is bigger or equal than the column index j.
		// The next statement fullfills this requirement.
		if(ea < eb) {
			std::swap(ea, eb);
		}
		entries.push_back(Entry{ea, eb, ew});
		diag[ea] += ew;
		diag[eb] += ew;
	}

	// do the conversion to a normal ev problem
	// assure global connectivity
	for(unsigned int i=0; i<diag.size(); i++) {
		Real& v = diag[i];
		if(v == 0) {
			// connect the disconnected cluster to all other clusters with a very small weight
			v = static_cast<Real>(1);
			Real q = static_cast<Real>(1) / static_cast<Real>(n-1);
			for(unsigned int j=0; j<i; j++) {
				auto it = std::find_if(entries.begin(), entries.end(), [i, j](const Entry& e) { return e.i == i && e.j == j; });
				if(it == entries.end()) {
					entries.push_back(Entry{i, j, q});
				}
			}
			for(unsigned int j=i+1; j<n; j++) {
				auto it = std::find_if(entries.begin(), entries.end(), [j, i](const Entry& e) { return e.i == j && e.j == i; });
				if(it == entries.end()) {
					entries.push_back(Entry{j, i, q});
				}
			}
			std::cerr << "Diagonal is 0! (i=" << i << ")" << std::endl;
		}
		else {
			v = static_cast<Real>(1) / std::sqrt(v);
		}
	}

	// a_ij for the transformed "normal" EV problem
	//		A x = \lambda x
	// is computed as follow from the diagonal matrix D and the weight
	// matrix W of the general EV problem
	//		(D - W) x = \lambda D x
	// as follows:
	//		a_ij = - w_ij / sqrt(d_i * d_j) if i != j
	//		a_ii = 1
	for(Entry& e : entries) {
		e.value = - e.value * diag[e.i] * diag[e.j];
	}
	for(unsigned int i=0; i<n; i++) {
		entries.push_back(Entry{i, i, static_cast<Real>(1)});
	}

	// sort entries to form a lower triangle matrix
	std::sort(entries.begin(), entries.end(), [](const Entry& a, const Entry& b) {
		return (a.j != b.j) ? (a.j < b.j) : (a.i < b.i);
	});



	typedef boost::numeric::ublas::symmetric_matrix<
		double, boost::numeric::ublas::lower> Matrix; 
	typedef boost::numeric::ublas::vector<double> Vector;

	int N = n;
	Matrix mat(N, N);
	for(int i=0;i<N;i++)
		for(int j=0;j<=i;j++)
			mat(i,j) = 0;   
	for(Entry& e : entries)
		mat(e.i,e.j) = e.value;

	typedef ietl::vectorspace<Vector> Vecspace;
	typedef boost::lagged_fibonacci607 Gen;  

	Vecspace vec(N);
	Gen mygen;
	ietl::lanczos<Matrix,Vecspace> lanczos(mat,vec);

	// Creation of an iteration object:  
	int max_iter = 10*N;  
	double rel_tol = 5000*std::numeric_limits<double>::epsilon();
	double abs_tol = 0.0001f;// std::pow(std::numeric_limits<double>::epsilon(),2./3);  
	std::cout << "Computation of 2 lowest converged eigenvalues\n\n";
	std::cout << "-----------------------------------\n\n";
	int n_lowest_eigenval = num_ev;
	std::vector<double> eigen;
	std::vector<double> err;
	std::vector<int> multiplicity;  
	ietl::lanczos_iteration_nlowest<double> iter(max_iter, n_lowest_eigenval, rel_tol, abs_tol);
	try{
		lanczos.calculate_eigenvalues(iter,mygen);
		//lanczos.more_eigenvalues(iter); 
		eigen = lanczos.eigenvalues();
		err = lanczos.errors();
		multiplicity = lanczos.multiplicities();
		std::cout<<"number of iterations: "<<iter.iterations()<<"\n";
	}
	catch (std::runtime_error& e) {
		std::cout << e.what() << "\n";
	} 

  // Printing eigenvalues with error & multiplicities:  
	std::cout << "#        eigenvalue            error         multiplicity\n";  
	std::cout.precision(10);
	for (int i=0;i<eigen.size();++i) 
		std::cout << i << "\t" << eigen[i] << "\t" << err[i] << "\t" << multiplicity[i] << "\n";

	// call of eigenvectors function follows:   
	std::cout << "\nEigen vectors computations for lowest eigenvalues:\n\n";  
	auto ew_begin = eigen.begin();
	while(*ew_begin <= 0.0f) ew_begin ++;
	auto ew_end = ew_begin + num_ev;
	std::vector<Vector> eigenvectors; // for storing the eigen vectors. 
	ietl::Info<double> info; // (m1, m2, ma, eigenvalue, residualm, status).

	try {
		lanczos.eigenvectors(ew_begin,ew_end,std::back_inserter(eigenvectors),info,mygen); 
	}
	catch (std::runtime_error& e) {
		std::cout << e.what() << "\n";
	}

	// std::cout << "Printing eigenvectors:\n\n"; 
	// for(std::vector<Vector>::iterator it = eigenvectors.begin();it!=eigenvectors.end();it++){
	// 	std::copy((it)->begin(),(it)->end(),std::ostream_iterator<double>(std::cout,", "));
	// 	std::cout << "\n\n";
	// }
	std::cout << " Information about the eigenvector computations:\n\n";
	for(int i = 0; i < info.size(); i++) {
		std::cout << " m1(" << i+1 << "): " << info.m1(i) << ", m2(" << i+1 << "): "
			<< info.m2(i) << ", ma(" << i+1 << "): " << info.ma(i) << " eigenvalue("
			<< i+1 << "): " << info.eigenvalue(i) << " residual(" << i+1 << "): "
			<< info.residual(i) << " error_info(" << i+1 << "): "
			<< info.error_info(i) << "\n\n";
	}


  	std::vector<EigenComponent> solution(eigenvectors.size());
	for(unsigned int i=0; i<solution.size(); i++) {
		EigenComponent& cmp = solution[i];
		cmp.eigenvalue = *(ew_begin + i);
#ifdef SPECTRAL_VERBOSE
		std::cout << "Eigenvalue " << i << ": " << cmp.eigenvalue << std::endl;
		cmp.eigenvector = Vec(n);
#endif
		for(unsigned int j=0; j<N; j++) {
			Real v = eigenvectors[i][j];
			// convert back to generalized eigenvalue problem!
			v *= diag[j];
			cmp.eigenvector[j] = v;
		}
	}
#ifdef SPECTRAL_VERBOSE
	std::cout << "Sparse Solver: returning" << std::endl;
#endif

	return solution;
#endif
}

}}

#endif
