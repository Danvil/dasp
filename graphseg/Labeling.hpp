/*
 * Segmentation.hpp
 *
 *  Created on: Mar 26, 2012
 *      Author: david
 */

#ifndef GRAPHSEG_LABELING_HPP_
#define GRAPHSEG_LABELING_HPP_

#include "Common.hpp"
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#include <vector>

namespace graphseg
{

	/** Node labeling using unique consecutive unsigned integers */
	struct GraphLabeling
	{
		// label for each node
		std::vector<unsigned int> labels;

		// number of unique labels
		unsigned int num_labels;

		/** Computes consecutive labels starting with 0 */
		void relabel();

		/** Computes labeling from non-consecutive labels */
		static GraphLabeling CreateClean(const std::vector<unsigned int>& labels);

	};

	/** Using connected components to compute graph segments */
	template<typename Graph>
	GraphLabeling ComputeSegmentLabels_ConnectedComponents(const Graph& graph, float threshold);

	/** Using UCM algorithm to compute graph segments */
	template<typename Graph>
	GraphLabeling ComputeSegmentLabels_UCM(const Graph& graph, float threshold);

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
	template<typename WeightedGraph>
	GraphLabeling ComputeSegmentLabels(const WeightedGraph& graph, float threshold, ComputeSegmentLabelsStrategy strategy=ComputeSegmentLabelsStrategies::ConnectedComponents) {
		switch(strategy) {
		default: case ComputeSegmentLabelsStrategies::ConnectedComponents:
			return ComputeSegmentLabels_ConnectedComponents(graph, threshold);
		case ComputeSegmentLabelsStrategies::UCM:
			return ComputeSegmentLabels_UCM(graph, threshold);
		}
	}

	template<typename WeightedGraph>
	GraphLabeling ComputeSegmentLabels_ConnectedComponents(const WeightedGraph& graph, float threshold)
	{
		// create a graph will all edges with costs >= threshold
		SpectralGraph cropped(boost::num_vertices(graph));
		for(auto eid : as_range(boost::edges(graph))) {
			float weight = boost::get(boost::edge_weight_t(), graph, eid);
			// take only edges with weight < threshold
			if(weight <= threshold) {
				boost::add_edge(boost::source(eid, graph), boost::target(eid, graph), cropped);
			}
		}
		// compute connected components
		std::vector<unsigned int> cluster_labels(boost::num_vertices(cropped));
		unsigned int num_labels = boost::connected_components(cropped, &cluster_labels[0]);
		// return result;
		return GraphLabeling{ cluster_labels, num_labels };
	}

	template<typename WeightedGraph>
	GraphLabeling ComputeSegmentLabels_UCM(const WeightedGraph& graph, float threshold)
	{
		// every superpixel is one region
		std::vector<unsigned int> cluster_labels(boost::num_vertices(graph));
		for(unsigned int i=0; i<cluster_labels.size(); i++) {
			cluster_labels[i] = i;
		}
		// get edge data
		struct Edge {
			unsigned int a, b;
			float weight;
		};
		std::vector<Edge> edges;
		for(auto eid : as_range(boost::edges(graph))) {
			edges.push_back(Edge{
				static_cast<unsigned int>(boost::source(eid, graph)),
				static_cast<unsigned int>(boost::target(eid, graph)),
				boost::get(boost::edge_weight_t(), graph, eid)
			});
		}
		// sort edges by weight
		std::sort(edges.begin(), edges.end(), [](const Edge& x, const Edge& y){ return x.weight < y.weight; });
		// cut one edge after another
		for(unsigned int k=0; k<edges.size(); k++) {
			if(edges[k].weight >= threshold) {
				break;
			}
			// join segments connected by edge
			unsigned int l_old = cluster_labels[edges[k].a];
			unsigned int l_new = cluster_labels[edges[k].b];
			std::replace(cluster_labels.begin(), cluster_labels.end(), l_old, l_new);
	#ifdef SEGS_DBG_SHOWGUI
			if(k % 10 == 0) {
				std::cout << "Removed edge with weight " << edges[k].cost << std::endl;
				// show ucm
				relabel();
				slimage::Image3ub vis = computeLabelImage(clusters);
				slimage::gui::Show("ucm", vis);
			}
	#endif
		}
		// create continuous labels
		return GraphLabeling::CreateClean(cluster_labels);
	}


}

#endif
