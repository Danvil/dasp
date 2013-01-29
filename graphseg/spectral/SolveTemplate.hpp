/*
 * SolveTemplate.hpp
 *
 *  Created on: Okt 20, 2012
 *      Author: david
 */

#ifndef DASP_SPECTRAL_SOLVETEMPLATE_HPP_
#define DASP_SPECTRAL_SOLVETEMPLATE_HPP_

#include "SolveDenseTemplate.hpp"
#include "Magma.hpp"
#include "SolveSparseTemplate.hpp"
#include "SolveLapackTemplate.hpp"
#include "SpectralIetl.hpp"

namespace graphseg
{
	namespace detail
	{
		/** Computes n smallest eigenvalues/-vectors for a graph */
		template<typename Graph>
		std::vector<EigenComponent> SolveTemplate(const Graph& graph, unsigned int num_ev, SpectralMethod method) {
			// pick one eigenvalue more because the first one is omitted
			switch(method) {
				default: case SpectralMethod::Eigen:
					return detail::SolveDenseTemplate(graph, num_ev + 1);
				case SpectralMethod::Magma:
					return detail::SolveSpectral_Magma(graph, num_ev + 1);
				case SpectralMethod::ArpackPPSparse:
					return detail::SolveSparseTemplate(graph, num_ev + 1);
				case SpectralMethod::Lapack:
					return detail::SolveLapackTemplate(graph, num_ev + 1);
				case SpectralMethod::Ietl:
					return detail::SpectralIetl(graph, num_ev + 1);
		}
		}
	}
}

#endif
