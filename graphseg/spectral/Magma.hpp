/*
 * Magma.hpp
 *
 *  Created on: Jan 28, 2012
 *      Author: david
 */

#ifndef GRAPHSEG_SPECTRAL_MAGMA_HPP_
#define GRAPHSEG_SPECTRAL_MAGMA_HPP_

#include "../Common.hpp"
#include "../as_range.hpp"
#include <boost/graph/adjacency_list.hpp>
#include <iostream>
#include <fstream>
#include <vector>

namespace graphseg { namespace detail {

std::vector<EigenComponent> MagmaImpl(const Eigen::MatrixXf& A, const Eigen::VectorXf& D, unsigned int num_ev);

template<typename Graph>
std::vector<EigenComponent> SolveSpectral_Magma(const Graph& graph, unsigned int num_ev)
{
	unsigned int dim = boost::num_vertices(graph);
#ifdef SPECTRAL_VERBOSE
	unsigned int num_edges = boost::num_edges(graph);
	std::cout << "SpectralSegmentation: dimension=" << dim
			<< ", num_edges=" << num_edges
			<< ", matrix non-zero elements = " << 100*static_cast<float>(2 * num_edges) / static_cast<float>(dim*dim) << "%" << std::endl;
	std::vector<int> nodes_with_no_connection;
#endif
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
#ifdef SPECTRAL_VERBOSE
			nodes_with_no_connection.push_back(i);
#endif
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
#ifdef SPECTRAL_VERBOSE
	if(!nodes_with_no_connection.empty()) {
		std::cout << "Nodes without connections (#=" << nodes_with_no_connection.size() << "): ";
		for(int i : nodes_with_no_connection) {
			std::cout << i << ", ";
		}
		std::cout << std::endl;
	}
#endif	
	// compute matrices D = diagonal(Di) and A = D - W
	Eigen::VectorXf D = Eigen::VectorXf::Zero(dim);
	Mat A = -W;
	for(unsigned int i=0; i<dim; i++) {
		Real di = Di[i];
		BOOST_ASSERT(di > static_cast<Real>(0));
		A(i,i) += di;
		D(i) = di;
	}

	return MagmaImpl(A, D, num_ev);
}

}}

#endif
