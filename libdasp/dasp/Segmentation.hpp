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
#include <boost/graph/adjacency_list.hpp>
#include <vector>

struct borderpixels_t {
  typedef boost::edge_property_tag kind;
};

namespace dasp
{

extern std::vector<slimage::Image3ub> cSegmentationDebug;

/** Computes a border image for a graph with edge weights and edge borderpixels
 * Graph requires edge_weight_t and borderpixels_t.
 */
template<typename Graph>
slimage::Image1f CreateBorderImage(unsigned int w, unsigned int h, const Graph& graph);

/** Computes a border image for a graph with edge weights and edge borderpixels
 * Graph requires edge_weight_t and borderpixels_t.
 */
template<typename Graph>
slimage::Image1f CreateBorderImageInv(unsigned int w, unsigned int h, const Graph& graph);

/** Smoothes a contour image and convertes to unsigned char */
slimage::Image1ub CreateSmoothedContourImage(const slimage::Image1f& src, float scl);

struct ClusterLabeling
{
	// one label for each superpixel
	std::vector<unsigned int> labels;

	// number of unique labels
	unsigned int num_labels;

	/** Computes consecutive labels starting with 0 */
	void relabel();

	static ClusterLabeling CreateClean(const std::vector<unsigned int>& labels);

};

/** Computes 0/1 segment boundaries from superpixel segment labels using edge detection */
slimage::Image1ub CreateBoundariesFromLabels(const Superpixels& clusters, const ClusterLabeling& cluster_labels);

/** Creates superpixel segment labels from boundaries with a weight bigger than a threshold using connected components
 * This functions uses connected components to perform a flood fill.
 * Graph requires edge_weight_t.
 * @return (labels, num_labels)
 */
template<typename Graph>
ClusterLabeling CreateLabelsFromBoundaries(const Superpixels& clusters, const Graph& graph, float threshold);

/** Creates a UCM labeling for a graph
 * Graph requires edge_weight_t.
 */
template<typename Graph>
ClusterLabeling UCM(const Superpixels& clusters, const Graph& graph, float threshold);

/** Creates colors for labeling */
std::vector<slimage::Pixel3ub> ComputeSegmentColors(const Superpixels& clusters, const ClusterLabeling& labeling);

/** Creates an image where each superpixel is colored with the corresponding label color */
slimage::Image3ub ComputeLabelImage(const Superpixels& clusters, const ClusterLabeling& labeling, const std::vector<slimage::Pixel3ub>& colors);

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

typedef boost::adjacency_list<boost::vecS,boost::vecS,boost::undirectedS,boost::no_property,boost::edge_weight_t> EdgeWeightGraph;

EdgeWeightGraph SpectralSegmentation(const Superpixels& clusters, const SpectralSettings& settings=SpectralSettings());

EdgeWeightGraph MinCutSegmentation(const Superpixels& clusters);

}

#include "SegmentationImpl.hpp"

#endif
