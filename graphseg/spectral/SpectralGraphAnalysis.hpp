/*
 * SpectralGraphAnalysis.hpp
 *
 *  Created on: Dec 30, 2012
 *      Author: david
 */

#ifndef DASP_SPECTRAL_SPECTRALGRAPHANALYSIS_HPP_
#define DASP_SPECTRAL_SPECTRALGRAPHANALYSIS_HPP_

#include "../Common.hpp"
#include "../as_range.hpp"
#include "SolveTemplate.hpp"
#include <boost/graph/adjacency_list.hpp>
#include <iostream>
#include <vector>
#ifdef SEGS_DBG_PRINT
#include <fstream>
#endif

namespace graphseg
{
	namespace detail
	{
		/** Assembles edge weights from eigenvalues and eigenvectors */
		template<typename Graph>
		Vec AssembleEdgeWeights(const Graph& graph, const std::vector<EigenComponent>& solution)
		{
			Vec edge_weight = Vec::Zero(boost::num_edges(graph));
		//	// later we weight by eigenvalues
		//	// find a positive eigenvalue (need to do this because of ugly instabilities ...
		//	Real ew_pos = -1.0f;
		//	for(unsigned int i=0; ; i++) {
		//		if(solver.eigenvalues()[i] > 0) {
		//			// F IXME magic to get a not too small eigenvalue
		////			unsigned int x = (n_used_ew + i)/2;
		//			unsigned int x = i + 5;
		//			ew_pos = solver.eigenvalues()[x];
		//			break;
		//		}
		//	}
		//	// compute normalized weights from eigenvalues
		//	Vec weights = Vec::Zero(n_used_ew);
		//	for(unsigned int k=0; k<n_used_ew; k++) {
		//		Real ew = solver.eigenvalues()[k + 1];
		//		if(ew <= ew_pos) {
		//			ew = ew_pos;
		//		}
		//		weights[k] = 1.0f / std::sqrt(ew);
		//	}
		//	std::cout << "Weights = " << weights.transpose() << std::endl;
			// look into first eigenvectors
			// skip first component
			for(unsigned int k=1; k<solution.size(); k++) {
				const EigenComponent& eigen = solution[k];
				// omit if eigenvalue is not positive
				Real ew = eigen.eigenvalue;
				// FIXME this is due to numerical instabilities
				if(ew <= Real(0)) {
					continue;
				}
				// weight by eigenvalue
				float w = 1.0f / std::sqrt(ew);
				// get eigenvector and normalize
				Vec ev = eigen.eigenvector;
				ev = (ev - ev.minCoeff() * Vec::Ones(ev.rows())) / (ev.maxCoeff() - ev.minCoeff());
				// for each edge compute difference of eigenvector values
				Vec e_k = Vec::Zero(edge_weight.rows());
				// FIXME proper edge indexing
				unsigned int eid_index = 0;
				for(auto eid : as_range(boost::edges(graph))) {
					e_k[eid_index] = std::abs(ev[boost::source(eid, graph)] - ev[boost::target(eid, graph)]);
					eid_index++;
				}
#ifdef SPECTRAL_VERBOSE
				std::cout << "w=" << w << " e_k.maxCoeff()=" << e_k.maxCoeff() << std::endl;
#endif
		//		e_k /= e_k.maxCoeff();
		//		for(unsigned int i=0; i<e_k.rows(); i++) {
		//			e_k[i] = std::exp(-e_k[i]);
		//		}
				e_k *= w;

#ifdef SEGS_DBG_PRINT
				{
					std::ofstream ofs((boost::format("/tmp/edge_weights_%03d.txt") % k).str());
					for(unsigned int i=0; i<e_k.rows(); i++) {
						ofs << e_k[i] << std::endl;
					}
				}
#endif
				//
				edge_weight += e_k;
			}
			return edge_weight;
		}
	}

	/** Applies graphseg graph theory fu to a weighted undirected graph */
	template<typename Graph>
	Graph SpectralGraphAnalysis(const Graph& graph, unsigned int num_ev, bool use_dense_solver) {
		// pick one more because the first one is omitted
		std::vector<detail::EigenComponent> solution = detail::SolveTemplate(graph, num_ev);
		detail::Vec weights = AssembleEdgeWeights(graph, solution);
		// create new graph with resulting edge strength
		Graph result = graph;
		// FIXME proper edge indexing
		unsigned int eid_index = 0;
		for(auto eid : as_range(boost::edges(result))) {
			boost::put(boost::edge_weight, result, eid, weights[eid_index]);
			eid_index++;
		}
		return result;
	}

}

#endif
