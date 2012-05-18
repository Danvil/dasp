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
#include <boost/graph/copy.hpp>
#include <fstream>
#include <iostream>
#include <set>

//#define SEGS_DBG_SHOWGUI
//#define SEGS_DBG_CREATE_EV_IMAGES
//#define SEGS_DBG_PRINT
//#define SEGS_VERBOSE

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
#ifdef SEGS_DBG_SHOWGUI
	{
		slimage::Image3ub vis = clusters.color_raw.clone();
		dasp::plots::PlotEdges(vis, clusters.ComputeLabels(), slimage::Pixel3ub{{255,255,255}},1);
		slimage::gui::Show("color", vis);
	}
	{
		slimage::Image3ub vis = dasp::plots::PlotClusters(clusters, dasp::plots::ClusterPoints, dasp::plots::Color);
		dasp::plots::PlotEdges(vis, clusters.ComputeLabels(), slimage::Pixel3ub{{255,255,255}},1);
		slimage::gui::Show("clusters", vis);
	}

#endif

	unsigned int num_vertices = boost::num_vertices(graph);
#ifdef SEGS_VERBOSE
	std::cout << "SpectralSegmentation: n = " << n << std::endl;
#endif

	using namespace detail;

	std::vector<Entry> entries;
	for(auto eid : as_range(boost::edges(graph))) {
		entries.push_back(Entry{
			source_superpixel_id(eid, graph),
			target_superpixel_id(eid, graph),
			boost::get(weights, eid)
		});
	}

#ifdef SEGS_VERBOSE
	std::cout << "Edve connectivity: min=" << edge_connectivity.minCoeff() << ", max=" << edge_connectivity.maxCoeff() << std::endl;
#endif

#ifdef SEGS_DBG_SHOWGUI
	{
		slimage::Image3ub edge_connectivity_img(clusters.width(), clusters.height(), slimage::Pixel3ub{{0,0,0}});
		for(unsigned int eid=0; eid<Gnb.edges.size(); eid++) {
			for(unsigned int pid : border_pixels[eid]) {
				edge_connectivity_img[pid] = dasp::plots::IntensityColor(static_cast<float>(edge_connectivity[eid]), 0.0f, 1.0f);
			}
		}
		slimage::gui::Show("edge_connectivity", edge_connectivity_img);
	}
#endif

	Vec result_ew;
	Mat result_ev;
	SolveSpectral(entries, num_vertices, result_ew, result_ev);


	unsigned int n_used_ew = std::min(num_vertices - 1, cNEV);
#ifdef SEGS_VERBOSE
	std::cout << "Eigenvalues = " << result_ew.topRows(n_used_ew + 1).transpose() << std::endl;
#endif
#ifdef SEGS_DBG_CREATE_EV_IMAGES
	{
		cSegmentationDebug.clear();
		// create image from eigenvectors (omit first)
		for(unsigned int k=0; k<std::min(cNEV,3u); k++) {
			// get k-th eigenvector
			Vec ev = result_ev.col(k + 1);
			// convert to plotable values
			std::vector<unsigned char> ev_ub(n);
			for(unsigned int i=0; i<n; i++) {
				float v = 0.5f + 2.0f*ev[i];
				ev_ub[i] = static_cast<unsigned char>(std::min(255, std::max(0, static_cast<int>(255.0f * v))));
			}
			// write to image
			slimage::Image3ub img(clusters.width(), clusters.height(), slimage::Pixel3ub{{255,0,0}});
			clusters.ForPixelClusters([&img,&ev_ub](unsigned int cid, const dasp::Cluster& c, unsigned int pid, const dasp::Point& p) {
				unsigned char v = ev_ub[cid];
				img[pid] = slimage::Pixel3ub{{v,v,v}};
			});
			cSegmentationDebug.push_back(img);
#ifdef SEGS_DBG_SHOWGUI
			slimage::gui::Show((boost::format("ev_%02d") % (k+1)).str(), img);
#endif
		}
	}	// DEBUG
#endif
	detail::Vec edge_weight = detail::Vec::Zero(num_vertices);
//	// later we weight by eigenvalues
//	// find a positive eigenvalue (need to do this because of ugly instabilities ...
//	Real ew_pos = -1.0f;
//	for(unsigned int i=0; ; i++) {
//		if(solver.eigenvalues()[i] > 0) {
//			// FIXME magic to get a not too small eigenvalue
////			unsigned int x = (n_used_ew + i)/2;
//			unsigned int x = i + 5;
//			ew_pos = solver.eigenvalues()[x];
//			break;
//		}
//	}
//	// compute normalized weights from eigenvalues
//	Vec weights = Vec::Zero(n_used_ew);
//	for(unsigned int k=0; k<n_used_ew; k++) {
//		Real ew = solver.eigenvalues()[k + 1];
//		if(ew <= ew_pos) {
//			ew = ew_pos;
//		}
//		weights[k] = 1.0f / std::sqrt(ew);
//	}
//	std::cout << "Weights = " << weights.transpose() << std::endl;
	// look into first eigenvectors
	for(unsigned int k=0; k<n_used_ew; k++) {
		// weight by eigenvalue
		Real ew = result_ew[k + 1];
		if(ew <= Real(0)) {
			// omit if eigenvalue is not positive
			continue;
		}
		float w = 1.0f / std::sqrt(ew);
		// get eigenvector and normalize
		Vec ev = result_ev.col(k + 1);
		ev = (ev - ev.minCoeff()*Vec::Ones(ev.rows())) / (ev.maxCoeff() - ev.minCoeff());
		// for each edge compute difference of eigenvector values
		Vec e_k = Vec::Zero(num_vertices);
		// FIXME proper edge indexing
		unsigned int eid_index = 0;
		for(auto eid : as_range(boost::edges(graph))) {
			e_k[eid_index] = std::abs(ev[source_superpixel_id(eid, graph)] - ev[target_superpixel_id(eid, graph)]);
			eid_index++;
		}
#ifdef SEGS_VERBOSE
		std::cout << "w=" << w << " e_k.maxCoeff()=" << e_k.maxCoeff() << std::endl;
#endif
//		e_k /= e_k.maxCoeff();
//		for(unsigned int i=0; i<e_k.rows(); i++) {
//			e_k[i] = std::exp(-e_k[i]);
//		}

		e_k *= w;

#ifdef SEGS_DBG_PRINT
		{
			std::ofstream ofs((boost::format("/tmp/edge_weights_%03d.txt") % k).str());
			for(unsigned int i=0; i<e_k.rows(); i++) {
				ofs << e_k[i] << std::endl;
			}
		}
#endif

		//
		edge_weight += e_k;
	}

#ifdef SEGS_DBG_PRINT
	{
		std::ofstream ofs("/tmp/edge_weights_sum.txt");
		for(unsigned int i=0; i<edge_weight.rows(); i++) {
			ofs << edge_weight[i] << std::endl;
		}
	}
#endif


#ifdef SEGS_VERBOSE
	std::cout << "Edge weights: min=" << edge_weight.minCoeff() << ", max=" << edge_weight.maxCoeff() << std::endl;
	std::cout << "Edge weights: median=" << edge_weight[edge_weight.rows()/2] << std::endl;
#endif

//	edge_weight /= edge_weight[(95*edge_weight.rows())/100];
//	edge_weight /= edge_weight.maxCoeff();
//	std::cout << "Edge weights = " << edge_weight.transpose() << std::endl;

//	// original edge connectivity graph
//	graph::Graph graph_original(Gnb.numNodes());
//	for(unsigned int eid=0; eid<Gnb.getEdges().size(); eid++) {
//		graph::Edge e = Gnb.getEdges()[eid];
//		e.cost = edge_connectivity[eid];
//		graph_original.add(e);
//	}
//
	// create superpixel neighbourhood graph with edge strength
	EdgeWeightGraph result = detail::CreateSuperpixelGraph<EdgeWeightGraph>(num_vertices);
	{
		// FIXME proper edge indexing
		unsigned int eid_index = 0;
		for(auto eid : as_range(boost::edges(graph))) {
			auto edge = boost::add_edge(source_superpixel_id(eid, graph), target_superpixel_id(eid, graph), result);
			boost::put(boost::edge_weight_t(), result, edge.first, edge_weight[eid_index]);
			eid_index++;
		}
	}

#ifdef SEGS_DBG_SHOWGUI
	slimage::Image1ub boundaries = CreateBorderImage(clusters.width(), clusters.height(), graph); // FIXME <- fuse border_pixels
	slimage::gui::Show("boundaries", segs.boundaries, 0.03f);
#endif

#ifdef SEGS_DBG_SHOWGUI
	slimage::gui::WaitForKeypress();
#endif

	return result;
}

}
