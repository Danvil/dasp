/*
 * Segmentation.cpp
 *
 *  Created on: Mar 26, 2012
 *      Author: david
 */

#include "Plots.hpp"
#include "graphseg/Spectral.hpp"
#include <Eigen/Dense>
//#include <Eigen/Sparse>
#include <boost/format.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/copy.hpp>
#include <fstream>
#include <iostream>
#include <set>

//#define SEGS_DBG_SHOWGUI
//#define SEGS_DBG_CREATE_EV_IMAGES
//#define SEGS_DBG_PRINT
//#define SEGS_VERBOSE

#ifdef SEGS_DBG_SHOWGUI
#	include <Slimage/Gui.hpp>
#endif
 
namespace dasp
{

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

template<typename SuperpixelGraph, typename WeightMap>
EdgeWeightGraph SpectralSegmentation(const SuperpixelGraph& graph, WeightMap weights, unsigned int cNEV)
{
	// create graph for spectral solving
	graphseg::SpectralGraph spectral;
	boost::copy_graph(graph, spectral,
			boost::edge_copy(
				[&spectral,&weights](typename SuperpixelGraph::edge_descriptor src, typename graphseg::SpectralGraph::edge_descriptor dst) {
					boost::put(boost::edge_weight, spectral, dst, boost::get(weights, src));
				}
	));
	// do spectral graph foo
	graphseg::SpectralGraph solved = graphseg::SolveSpectral(spectral, cNEV);
	// create superpixel neighbourhood graph with edge strength
	EdgeWeightGraph result;
	boost::copy_graph(solved, result,
			boost::edge_copy(
				[&solved,&result](typename graphseg::SpectralGraph::edge_descriptor src, typename EdgeWeightGraph::edge_descriptor dst) {
					boost::put(boost::edge_weight, result, dst, boost::get(boost::edge_weight, solved, src));
				}
	));
	return result;
}

}
