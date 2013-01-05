/*
 * Segmentation.hpp
 *
 *  Created on: Mar 26, 2012
 *      Author: david
 */

#ifndef GRAPHSEG_LABELING_HPP_
#define GRAPHSEG_LABELING_HPP_

#include "Common.hpp"
#include "as_range.hpp"
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
		// cut one edge after another and merge segments
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

	template<typename Graph, typename EdgeWeightMap, typename VertexLabelMap>
	inline GraphLabeling ComputeSegmentLabels_UCM_Supervised(
		const Graph& graph,
		EdgeWeightMap edge_weight_map,
		VertexLabelMap vertex_label_map,
		float threshold)
	{
		const std::size_t num_vertices = boost::num_vertices(graph);
		// find maximal used label
		int label_max = -1;
		for(auto vid : as_range(boost::vertices(graph))) {
			const int label = boost::get(vertex_label_map, vid);
			label_max = std::max(label_max, label);
		}
		const int label_supervised_threshold = label_max;
		// labeled clusters keep label, unlabeled clusters get new label
		// labels will indicate if clusters are supervised:
		//   is_supervised(c) := (label(c) <= label_supervised_threshold)
		std::vector<unsigned int> cluster_labels(num_vertices);
		for(auto vid : as_range(boost::vertices(graph))) {
			const int label = graph[vid];
			if(label != -1) {
				// supervised cluster
				cluster_labels[vid] = label;
			}
			else {
				// unsupervised cluster
				cluster_labels[vid] = ++label_max;
			}
		}
		// get all edges
		struct Edge
		{
			unsigned int a, b;
			float weight;
		};
		std::vector<Edge> edges;
		edges.reserve(boost::num_edges(graph));
		for(auto eid : as_range(boost::edges(graph))) {
			edges.push_back(Edge{
				static_cast<unsigned int>(boost::source(eid, graph)),
				static_cast<unsigned int>(boost::target(eid, graph)),
				boost::get(edge_weight_map, eid)
			});
		}
		// sort edges by weight
		std::sort(edges.begin(), edges.end(), [](const Edge& x, const Edge& y) { return x.weight < y.weight; });
		// cut one edge after another and merge clusters into segments
		for(unsigned int k=0; k<edges.size(); k++) {
			const Edge& edge = edges[k];
			// merge only if edge weight is smaller than thresold
			if(edge.weight >= threshold) {
				// all consecutive edges will have larger edge weights -> break
				break;
			}
			// get cluster labels
			const unsigned int label_a = cluster_labels[edge.a];
			const unsigned int label_b = cluster_labels[edge.b];
			if(label_a == label_b) {
				// clusters are already part of the same segment -> nothing to do
				continue;
			}
			// check if cluster segments are supervised
			const bool is_supervised_a = (label_a <= label_supervised_threshold);
			const bool is_supervised_b = (label_b <= label_supervised_threshold);
			if(is_supervised_a && is_supervised_b) {
				// do not merge if both clusters are supervised
				continue;
			}
			// join segments connected by edge and keep the smaller label
			// this assures that the label supervision property is guaranteed
			const unsigned int old_label = std::max(label_a, label_b);
			const unsigned int new_label = std::min(label_a, label_b);
			std::replace(cluster_labels.begin(), cluster_labels.end(), old_label, new_label);
		}
		// create continuous labels
		return GraphLabeling::CreateClean(cluster_labels);
	}

}

#endif
