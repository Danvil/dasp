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

	/** Computes 0/1 segment boundaries from superpixel segment labels using edge detection */
	void createBoundariesFromLabels(const Clustering& clusters);

	/** Creates superpixel segment labels from boundaries with a weight bigger than a threshold using connected components */
	void createLabelsFromBoundaries(const Clustering& clusters, float threshold);

	std::vector<slimage::Pixel3ub> computeSegmentColors(const Clustering& clusters) const;

	slimage::Image3ub computeLabelImage(const Clustering& clusters) const;

};

Segmentation SpectralSegmentation(const Clustering& clusters);

Segmentation MinCutSegmentation(const Clustering& clusters);

}

#endif
