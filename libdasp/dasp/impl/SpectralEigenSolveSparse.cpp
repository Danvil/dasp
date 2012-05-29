/*
 * SpectralEigenSolveSparse.cpp
 *
 *  Created on: May 29, 2012
 *      Author: david
 */

#include "Spectral.hpp"
#include <arpack++/arlssym.h>
#include <cmath>
#include <iostream>

void MemoryOverflow() {
	std::cerr << "ArpackError: MEMORY_OVERFLOW" << std::endl;
	throw ArpackError(ArpackError::MEMORY_OVERFLOW);
}

void ArpackError::Set(ArpackError::ErrorCode code, char const* msg) {
//	ArpackError::code = code;
	std::cerr << "ArpackError: code=" << code << ", msg=" << std::string(msg) << std::endl;
}

namespace dasp { namespace detail {

PartialEigenSolution SpectralEigenSolveSparse(const SpectralGraph& graph, unsigned int num_ev)
{
	std::cout << "Sparse Solver: started" << std::endl;

	// We want to solve the EV problem: (D - W) x = \lamda D x.
	// Each edge of the graph defines two entries into the symmetric matrix W.
	// The diagonal matrix D is defined via d_i = sum_j{w_ij}.

	// As D is a diagonal matrix the the general problem can be easily transformed
	// into a normal eigenvalue problem by decomposing D = L L^t, which yields L = sqrt(D).
	// Thus the EV problem is: L^{-1} (D - W) L^{-T} y = \lambda y.
	// Eigenvectors can be transformed using x = L^{-T} y.

	// The dimension of the problem
	int n = boost::num_vertices(graph);

	// Each edge defines two entries (one in the upper and one in the lower).
	// In addition all diagonal entries are non-zero.
	// Thus the number of non-zero entries in the lower triangle is equal to
	// the number of edges plus the number of nodes.
	int nnz = boost::num_edges(graph) + n;
	// This vector will hold the non-zero elements of the lower triangle.
	std::vector<Real> nzval(nnz);

	// collect all non-zero elements
	struct Entry {
		int i, j;
		Real value;
	};
	std::vector<Entry> entries;
	entries.reserve(nnz);

	// also collect diagonal entries
	std::vector<Real> diag(n);

	// no collect entries
	for(auto eid : as_range(boost::edges(graph))) {
		int ea = static_cast<int>(boost::source(eid, graph));
		int eb = static_cast<int>(boost::target(eid, graph));
		float ew = boost::get(boost::edge_weight, graph, eid);
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
	for(Real& v : diag) {
		v = static_cast<Real>(1) / std::sqrt(v);
	}
	for(Entry& e : entries) {
		e.value *= (diag[e.i] * diag[e.j]);
	}

	// add the diagonal entries which are are all 1 (for the normal EV problem)
	for(unsigned int i=0; i<n; i++) {
		entries.push_back(Entry{i, i, static_cast<Real>(1)});
	}

	// sort entries to form a lower triangle matrix
	std::sort(entries.begin(), entries.end(), [](const Entry& a, const Entry& b) {
		return (a.j != b.j) ? (a.j < b.j) : (a.i < b.i);
	});

	// read out structure
	// Assumes that each column has at least one non-zero element.
	std::vector<int> irow(nnz);
	std::vector<int> pcol;
	pcol.reserve(n + 1);
	int current_col = -1;
	for(unsigned int i=0; i<entries.size(); i++) {
		const Entry& e = entries[i];
		irow[i] = e.i;
		if(e.j == current_col + 1) {
			pcol.push_back(i);
			current_col++;
		}
	}
	pcol.push_back(nnz);

	// define ARPACK matrix
	std::cout << "Sparse Solver: defining matrix" << std::endl;
	ARluSymMatrix<Real> mat(n, nnz, nzval.data(), irow.data(), pcol.data());

	// solve ARPACK problem
	std::cout << "Sparse Solver: solving ..." << std::flush;
	ARluSymStdEig<Real> solv(num_ev, mat);
	std::vector<Real> v_ew(num_ev);
	std::vector<Real> v_ev(num_ev * n);
	Real* p_ew = v_ew.data();
	Real* p_ev = v_ev.data();
	solv.EigenValVectors(p_ev, p_ew, false);
	std::cout << " finished." << std::endl;

	std::cout << "Sparse Solver: collecting results" << std::endl;

	PartialEigenSolution solution(num_ev);
	for(unsigned int i=0; i<num_ev; i++) {
		EigenComponent& cmp = solution[i];
		cmp.eigenvalue = p_ew[i];
		std::cout << "Eigenvalue " << i << ": " << cmp.eigenvalue << std::endl;
		cmp.eigenvector = Vec(n);
		for(unsigned int j=0; j<n; j++) {
			cmp.eigenvector[j] = p_ev[i*n + j];
		}
		// FIXME convert back !!!
	}

	std::cout << "Sparse Solver: returning" << std::endl;
	return solution;

	//return SpectralDenseEigenSolve(graph, num_ev);

	// FIXME

	//	{	// DEBUG
	//		std::ofstream ofs_D("/tmp/spectral_D.csv");
	//		std::ofstream ofs_W("/tmp/spectral_W.csv");
	//		for(unsigned int i=0; i<n; i++) {
	//			for(unsigned int j=0; j<n; j++) {
	//				ofs_D << D(i,j);
	//				ofs_W << W(i,j);
	//				if(j+1 == n) {
	//					ofs_D << std::endl;
	//					ofs_W << std::endl;
	//				}
	//				else {
	//					ofs_D << ",";
	//					ofs_W << ",";
	//				}
	//			}
	//		}
	//	}	// DEBUG
}

}}
