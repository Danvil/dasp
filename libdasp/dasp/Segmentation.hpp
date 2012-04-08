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

	graph::Graph original_graph;

	graph::Graph segmentation_graph;

	slimage::Image1f original_boundaries;

	slimage::Image1f boundaries;

	slimage::Image1ub boundaries_wt;

	unsigned int countSegments() const {
		return segment_count;
	}

	/** Computes 0/1 segment boundaries from superpixel segment labels using edge detection */
	void createBoundariesFromLabels(const Superpixels& clusters);

	/** Creates superpixel segment labels from boundaries with a weight bigger than a threshold using connected components */
	void createLabelsFromBoundaries(const Superpixels& clusters, float threshold);

	void ucm(const Superpixels& clusters, float threshold);

	std::vector<slimage::Pixel3ub> computeSegmentColors(const Superpixels& clusters) const;

	slimage::Image3ub computeLabelImage(const Superpixels& clusters) const;

	void relabel();

	static slimage::Image1f CreateBorderImage(unsigned int w, unsigned int h, const graph::Graph& graph, const std::vector<std::vector<unsigned int>>& border_pixels);

	static slimage::Image1f CreateBorderImageInv(unsigned int w, unsigned int h, const graph::Graph& graph, const std::vector<std::vector<unsigned int>>& border_pixels);

	static slimage::Image1ub CreateSmoothedContourImage(const slimage::Image1f& src, float scl);

};

struct SpectralSettings
{
	SpectralSettings() {
		num_eigenvectors = 16;
		w_spatial = 1;
		w_normal = 1;
		w_color = 1;
	}
	unsigned int num_eigenvectors;
	float w_spatial;
	float w_normal;
	float w_color;
};

Segmentation SpectralSegmentation(const Superpixels& clusters, const SpectralSettings& settings=SpectralSettings());

Segmentation MinCutSegmentation(const Superpixels& clusters);

}

#endif
