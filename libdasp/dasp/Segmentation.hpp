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

/** Computes a border image
 * All border pixels for each edge are set to the edge weight in the image.
 * Default value for non border pixels is 0.
 * Border pixels are given in index form x + y*w where w must be the same given as parameter.
 * @param w width of result image
 * @param h width of result image
 * @param graph a boost graph
 * @param weights an weight property map of type float
 * @param border_pixels an edge property map of type std::vector<unsigned int>
 * @param return image with painted border pixels
 */
template<typename Graph, typename WeightMap, typename BorderPixelMap>
slimage::Image1f CreateBorderPixelImage(unsigned int w, unsigned int h, const Graph& graph, WeightMap weights, BorderPixelMap border_pixels);

/** Smoothes a contour image and converts to unsigned char */
slimage::Image1ub CreateSmoothedContourImage(const slimage::Image1f& src, float scl);

/** Node labeling using unique consecutive unsigned integers */
struct ClusterLabeling
{
	// label for each node
	std::vector<unsigned int> labels;

	// number of unique labels
	unsigned int num_labels;

	/** Computes consecutive labels starting with 0 */
	void relabel();

	/** Computes labeling from non-consecutive labels */
	static ClusterLabeling CreateClean(const std::vector<unsigned int>& labels);

};

namespace impl
{
	/** Using connected components to compute graph segments */
	template<typename Graph>
	ClusterLabeling ComputeSegmentLabels_ConnectedComponents(const Graph& graph, float threshold);

	/** Using UCM algorithm to compute graph segments */
	template<typename Graph>
	ClusterLabeling ComputeSegmentLabels_UCM(const Graph& graph, float threshold);
}

namespace ComputeSegmentLabelsStrategies
{
	enum type {
		ConnectedComponents,
		UCM
	};
}
typedef ComputeSegmentLabelsStrategies::type ComputeSegmentLabelsStrategy;

/** Computes graph segments separated by edges with big weight
 * Connects all vertices which are connected by edges with a weight smaller than the threshold.
 * Graph requires edge_weight_t.
 */
template<typename Graph>
ClusterLabeling ComputeSegmentLabels(const Graph& graph, float threshold, ComputeSegmentLabelsStrategy strategy=ComputeSegmentLabelsStrategies::ConnectedComponents) {
	switch(strategy) {
	default: case ComputeSegmentLabelsStrategies::ConnectedComponents: return impl::ComputeSegmentLabels_ConnectedComponents(graph, threshold);
	case ComputeSegmentLabelsStrategies::UCM: return impl::ComputeSegmentLabels_UCM(graph, threshold);
	}
}

/** Computes 0/1 segment boundaries from superpixel segment labels using edge detection */
slimage::Image1ub CreateBoundaryImageFromLabels(const Superpixels& clusters, const ClusterLabeling& cluster_labels);

/** Creates colors for labeling */
std::vector<slimage::Pixel3ub> ComputeSegmentColors(const Superpixels& clusters, const ClusterLabeling& labeling);

/** Creates an image where each superpixel is colored with the corresponding label color */
slimage::Image3ub CreateLabelImage(const Superpixels& clusters, const ClusterLabeling& labeling, const std::vector<slimage::Pixel3ub>& colors);

/** Graph type of segmentation result */
typedef boost::adjacency_list<boost::vecS,boost::vecS,boost::undirectedS,boost::no_property,boost::property<boost::edge_weight_t,float>> EdgeWeightGraph;

/** Settings for spectral graph segmentation */
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

/** Performs spectral graph segmentation */
EdgeWeightGraph SpectralSegmentation(const Superpixels& clusters, const SpectralSettings& settings=SpectralSettings());

/** Performs min-cut graph segmentation */
EdgeWeightGraph MinCutSegmentation(const Superpixels& clusters);

}

#include "SegmentationImpl.hpp"

#endif
