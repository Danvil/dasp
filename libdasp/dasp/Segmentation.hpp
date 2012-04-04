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
	void createBoundariesFromLabels(const Superpixels& clusters);

	/** Creates superpixel segment labels from boundaries with a weight bigger than a threshold using connected components */
	void createLabelsFromBoundaries(const Superpixels& clusters, float threshold);

	std::vector<slimage::Pixel3ub> computeSegmentColors(const Superpixels& clusters) const;

	slimage::Image3ub computeLabelImage(const Superpixels& clusters) const;

};

Segmentation SpectralSegmentation(const Superpixels& clusters);

Segmentation MinCutSegmentation(const Superpixels& clusters);

}

#endif
