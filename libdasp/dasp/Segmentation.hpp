/*
 * Segmentation.hpp
 *
 *  Created on: Mar 26, 2012
 *      Author: david
 */

#ifndef SEGMENTATION_HPP_
#define SEGMENTATION_HPP_

#include <Slimage/Slimage.hpp>
#include "Superpixels.hpp"
#include <vector>

namespace dasp
{

extern std::vector<slimage::Image3ub> cSegmentationDebug;

struct Segmentation
{
	std::vector<unsigned int> cluster_labels;

	unsigned int segment_count;

	graph::Graph segmentation_graph;

	slimage::Image1f boundaries;

	slimage::Image1ub boundaries_wt;

	unsigned int countSegments() const {
		return segment_count;
	}

	void createBoundariesFromLabels(const Clustering& clusters);

};

Segmentation SpectralSegmentation(const Clustering& clusters);

Segmentation MinCutSegmentation(const Clustering& clusters);

}

#endif
