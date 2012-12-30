/*
 * Graph.hpp
 *
 *  Created on: May 18, 2012
 *      Author: david
 */

#ifndef DASP_GRAPH_HPP_
#define DASP_GRAPH_HPP_

#include "graphseg/as_range.hpp"
#include <boost/graph/adjacency_list.hpp>
#include <Eigen/Dense>

namespace dasp
{
	struct borderpixels_t {
		typedef boost::edge_property_tag kind;
	};

	typedef boost::property<boost::edge_weight_t, float> EdgeWeightProperty;

	typedef boost::property<borderpixels_t, std::vector<unsigned int>> EdgeBorderPixelsProperty;

	/** An undirected graph with a list of image border pixels per edge */
	typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS,
		boost::no_property,
		EdgeBorderPixelsProperty
	> BorderPixelGraph;

	/** An undirected weighted graph */
	typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS,
		boost::no_property,
		EdgeWeightProperty
	> EdgeWeightGraph;

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
		EdgeWeightProperty
	> DaspGraph;

}

#endif
