/*
 * SuperpixelGraph.cpp
 *
 *  Created on: Mar 7, 2012
 *      Author: david
 */

#include "SuperpixelGraph.hpp"
#include "TreeReduction.hpp"
//----------------------------------------------------------------------------//
namespace dasp {
//----------------------------------------------------------------------------//

void SuperpixelGraph::createConnections(float threshold)
{
	std::size_t n = size();
	node_connections_.resize(n);
	for(std::size_t i=0; i<n; i++) {
		for(std::size_t j=i+1; j<n; j++) {
			float d = (nodes_[i].position - nodes_[j].position).norm();
			// only connect if distance is smaller than threshold
			if(d < threshold) {
				node_connections_[i].push_back(j);
				node_connections_[j].push_back(i);
			}
		}
	}
}

SuperpixelNeighbourhoodGraph SuperpixelGraph::createNeighbourhoodGraph(unsigned int i) const
{
	SuperpixelNeighbourhoodGraph ng;
	ng.center_ = nodes_[i];
	for(std::size_t j : node_connections_[i]) {
		ng.neighbours_.push_back(nodes_[j]);
	}
	return ng;
}

//----------------------------------------------------------------------------//
}
//----------------------------------------------------------------------------//
