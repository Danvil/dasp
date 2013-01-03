/*
 * Graph.hpp
 *
 *  Created on: May 18, 2012
 *      Author: david
 */

#ifndef DASP_GRAPH_HPP_
#define DASP_GRAPH_HPP_

#include "graphseg/as_range.hpp"
#include "graphseg/Spectral.hpp"
#include <boost/graph/adjacency_list.hpp>
#include <Eigen/Dense>

namespace dasp
{
	typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> UndirectedGraph;

	typedef graphseg::SpectralGraph UndirectedWeightedGraph;

	/** Core superpoint information */
	struct DaspPoint
	{
		Eigen::Vector2f px;
		Eigen::Vector3f position;
		Eigen::Vector3f color;
		Eigen::Vector3f normal;
	};

	/** Weighted graph of superpoints */
	typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS,
		DaspPoint,
		boost::property<boost::edge_weight_t, float>
	> DaspGraph;

}

#endif
