/*
 * SuperpixelGraph.hpp
 *
 *  Created on: Mar 7, 2012
 *      Author: david
 */

#ifndef SUPERPIXELGRAPH_HPP_
#define SUPERPIXELGRAPH_HPP_
//----------------------------------------------------------------------------//
#include <eigen3/Eigen/Dense>
#include <vector>
//----------------------------------------------------------------------------//
namespace dasp {
//----------------------------------------------------------------------------//

struct SuperpixelState
{
	unsigned int x, y;
	Eigen::Vector3f position;
	Eigen::Vector3f color;
	Eigen::Vector3f normal;
	float scala;
};

struct SuperpixelNeighbourhoodGraph
{
	SuperpixelState center_;
	std::vector<SuperpixelState> neighbours_;
};

struct SuperpixelGraph
{
	std::vector<SuperpixelState> nodes_;

	std::vector<std::vector<std::size_t>> node_connections_;

	std::size_t size() const {
		return nodes_.size();
	}

	void createConnections(float threshold);

	SuperpixelNeighbourhoodGraph createNeighbourhoodGraph(unsigned int i) const;

};

//----------------------------------------------------------------------------//
}
//----------------------------------------------------------------------------//
#endif
