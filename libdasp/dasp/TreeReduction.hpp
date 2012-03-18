/*
 * TreeReduction.hpp
 *
 *  Created on: Jan 29, 2012
 *      Author: david
 */

#ifndef TREEREDUCTION_HPP_
#define TREEREDUCTION_HPP_

#include <vector>

namespace Romeo {

namespace TreeReduction
{

	struct Edge {
		unsigned int a, b;
		float cost;
	};

	struct Graph {
		Graph() : nodes_(0) {}
		Graph(unsigned int n) : nodes_(n) {}
		unsigned int nodes_;
		std::vector<Edge> edges;
	};

	Graph MinimalCostEdges(const Graph& input);

	Graph MinimalSpanningCutting(const Graph& input, const float cut_param, std::vector<unsigned int>* labels=0);

};

}

#endif /* TREEREDUCTION_HPP_ */
