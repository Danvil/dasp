/*
 * SolveTemplate.hpp
 *
 *  Created on: Okt 20, 2012
 *      Author: david
 */

#ifndef DASP_SPECTRAL_SOLVETEMPLATE_HPP_
#define DASP_SPECTRAL_SOLVETEMPLATE_HPP_

#include "SolveDenseTemplate.hpp"
#include "SolveSparseTemplate.hpp"

namespace graphseg
{
	namespace detail
	{
		/** Computes n smallest eigenvalues/-vectors for a graph */
		template<typename Graph>
		std::vector<EigenComponent> SolveTemplate(const Graph& graph, unsigned int num_ev, bool use_dense_solver=true) {
			// pick one eigenvalue more because the first one is omitted
			if(use_dense_solver) {
				return detail::SolveDenseTemplate(graph, num_ev + 1);
			}
			else {
				return detail::SolveSparseTemplate(graph, num_ev + 1);
			}
		}
	}
}

#endif
