/*
 * Segmentation.hpp
 *
 *  Created on: Mar 26, 2012
 *      Author: david
 */

#ifndef SEGMENTATION_HPP_
#define SEGMENTATION_HPP_

#include "Superpixels.hpp"
#include "Graph.hpp"
#include <Slimage/Slimage.hpp>
#include <Eigen/Dense>
#include <boost/graph/adjacency_list.hpp>
#include <vector>

namespace dasp
{

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

/** Performs spectral graph segmentation */
template<typename SuperpixelGraph, typename WeightMap>
UndirectedWeightedGraph SpectralSegmentation(const SuperpixelGraph& graph, WeightMap weights, unsigned int num_eigenvectors=24);

}

#include "SegmentationImpl.hpp"

#endif
