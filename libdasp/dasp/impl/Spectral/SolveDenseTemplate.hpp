/*
 * SolveDenseTemplate.hpp
 *
 *  Created on: Okt 20, 2012
 *      Author: david
 */

#ifndef DASP_IMPL_SPECTRAL_SOLVEDENSETEMPLATE_HPP_
#define DASP_IMPL_SPECTRAL_SOLVEDENSETEMPLATE_HPP_

#include "Types.hpp"
#include <iostream>

namespace dasp { namespace detail {

template<typename Graph>
PartialEigenSolution SolveDenseTemplate(const Graph& graph, unsigned int num_ev)
{
	unsigned int dim = boost::num_vertices(graph);
	if(cVerbose) {
		unsigned int num_edges = boost::num_edges(graph);
		std::cout << "SpectralSegmentation: dimension=" << dim
				<< ", num_edges=" << num_edges
				<< ", matrix non-zero elements = " << 100*static_cast<float>(2 * num_edges) / static_cast<float>(dim*dim) << "%" << std::endl;
	}
	// creating matrices
	Mat W = Mat::Zero(dim,dim);
	std::vector<float> Di(dim, 0.0f);
	for(auto eid : as_range(boost::edges(graph))) {
		unsigned int ea = boost::source(eid, graph);
		unsigned int eb = boost::target(eid, graph);
		float ew = boost::get(boost::edge_weight, graph, eid);
		if(std::isnan(ew)) {
			std::cerr << "Weight for edge (" << ea << "," << eb << ") is nan!" << std::endl;
			continue;
		}
		if(ew < 0) {
			std::cerr << "Weight for edge (" << ea << "," << eb << ") is negative!" << std::endl;
			continue;
		}
		W(ea, eb) = ew;
		W(eb, ea) = ew;
		Di[ea] += ew;
		Di[eb] += ew;
	}
	// connect disconnected segments to everything
	// FIXME why is this necessary?
	for(unsigned int i=0; i<dim; i++) {
		float& di = Di[i];
		if(di == 0) {
			if(cVerbose)
				std::cout << "Node " << i << " has no connections! " << std::endl;
			// connect the disconnected cluster to all other clusters with a very small weight
			di = 1.0f;
			float q = di / static_cast<float>(dim-1);
			for(unsigned int j=0; j<dim; j++) {
				if(j == i) continue;
				W(i,j) = q;
				W(j,i) = q;
			}
		}
	}
	// compute matrices D = diagonal(Di) and A = D - W
	Mat D = Mat::Zero(dim,dim);
	Mat A = -W;
	for(unsigned int i=0; i<dim; i++) {
		Real di = Di[i];
		BOOST_ASSERT(di > static_cast<Real>(0));
		A(i,i) += di;
		D(i,i) = di;
	}
	// solve eigensystem
	Eigen::GeneralizedSelfAdjointEigenSolver<Mat> solver;
	solver.compute(A, D);
	if(cVerbose)
		std::cout << "SpectralSegmentation: GeneralizedSelfAdjointEigenSolver says " << solver.info() << std::endl;
	// return eigenvectors and eigenvalues
	PartialEigenSolution solution(std::min(dim, num_ev));
	for(unsigned int i=0; i<solution.size(); i++) {
		solution[i].eigenvalue = solver.eigenvalues()[i];
		if(cVerbose)
			std::cout << "SpectralSegmentation: eigenvalue #" << i << "=" << solution[i].eigenvalue << std::endl;
		solution[i].eigenvector = solver.eigenvectors().col(i);
	}
	return solution;
}

}}

#endif
