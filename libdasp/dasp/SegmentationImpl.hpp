/*
 * Segmentation.cpp
 *
 *  Created on: Mar 26, 2012
 *      Author: david
 */

#include "Plots.hpp"
#include <Slimage/Gui.hpp>
#include <Eigen/Dense>
//#include <Eigen/Sparse>
#include <Eigen/Eigenvalues>
#include <boost/assert.hpp>
#include <boost/format.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#include <fstream>
#include <iostream>
#include <set>

//#define SEGS_DBG_SHOWGUI
//#define SEGS_DBG_CREATE_EV_IMAGES
//#define SEGS_DBG_PRINT
//#define SEGS_VERBOSE

namespace dasp
{

template<typename Graph, typename WeightMap, typename BorderPixelMap>
slimage::Image1f CreateBorderPixelImage(unsigned int w, unsigned int h, const Graph& graph, WeightMap weights, BorderPixelMap border_pixels)
{
	slimage::Image1f result(w, h, slimage::Pixel1f{0.0f});
	for(auto eid : as_range(boost::edges(graph))) {
		float v = boost::get(weights, eid);
		for(unsigned int pid : boost::get(border_pixels, eid)) {
			result[pid] = v;
		}
	}
	return result;
}

template<typename Graph>
ClusterLabeling impl::ComputeSegmentLabels_ConnectedComponents(const Graph& graph, float threshold)
{
	// create a graph will all edges with costs >= threshold
	typedef boost::adjacency_list<boost::vecS,boost::vecS,boost::undirectedS,boost::no_property, boost::edge_weight_t> graph_t;
	graph_t cropped;
	for(unsigned int i=0; i<boost::num_vertices(graph); i++) {
		boost::add_vertex(cropped);
	}
	for(auto eid : as_range(boost::edges(graph))) {
		float weight = boost::get(boost::edge_weight_t(), graph, eid);
		// take only edges with weight < threshold
		if(weight <= threshold) {
			boost::add_edge(
					static_cast<unsigned int>(boost::source(eid, graph)),
					static_cast<unsigned int>(boost::target(eid, graph)),
					cropped);
		}
	}
	// compute connected components
	std::vector<unsigned int> cluster_labels(boost::num_vertices(cropped));
	unsigned int num_labels = boost::connected_components(cropped, &cluster_labels[0]);
	// return result;
	return ClusterLabeling{ cluster_labels, num_labels };
}

template<typename Graph>
ClusterLabeling impl::ComputeSegmentLabels_UCM(const Graph& graph, float threshold)
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
	return ClusterLabeling::CreateClean(cluster_labels);
}

}
