/*
 * Segmentation.cpp
 *
 *  Created on: Mar 26, 2012
 *      Author: david
 */

#include "Segmentation.hpp"
#include "tools/Graph.hpp"
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

std::vector<slimage::Image3ub> cSegmentationDebug;

slimage::Image1f Segmentation::CreateBorderImage(unsigned int w, unsigned int h, const graph::Graph& graph, const std::vector<std::vector<unsigned int>>& border_pixels)
{
	slimage::Image1f result(w, h, slimage::Pixel1f{0.0f});
	graph.foreachEdge([&border_pixels,&result](const graph::Edge& e) {
		float v = e.cost;
		for(unsigned int pid : border_pixels[e.id]) {
			result[pid] = v;
		}
	});
	return result;
}

slimage::Image1f Segmentation::CreateBorderImageInv(unsigned int w, unsigned int h, const graph::Graph& graph, const std::vector<std::vector<unsigned int>>& border_pixels)
{
	slimage::Image1f result(w, h, slimage::Pixel1f{0.0f});
	graph.foreachEdge([&border_pixels,&result](const graph::Edge& e) {
		float v = e.cost;
		for(unsigned int pid : border_pixels[e.id]) {
			result[pid] = 1.0f - v;
		}
	});
	return result;
}

slimage::Image1ub Segmentation::CreateSmoothedContourImage(const slimage::Image1f& src, float scl)
{
	slimage::Image1f tmp(src.dimensions(), slimage::Pixel1f{1.0f});
	for(unsigned int x=1; x<tmp.width()-1; x++) {
		for(unsigned int y=1; y<tmp.height()-1; y++) {
			float nb = src(x-1,y) + src(x+1,y) + src(x,y-1) + src(x,y+1);
			float v = src(x,y) + nb / 4.0f;
			tmp(x,y) = 1.0f - scl * v;
		}
	}
	slimage::Image1ub vis(tmp.dimensions());
	slimage::conversion::Convert(tmp, vis);
	return vis;
}

void Segmentation::createBoundariesFromLabels(const Superpixels& clustering)
{
	// create segment labeling
	slimage::Image1i labels(clustering.width(), clustering.height(), slimage::Pixel1i{-1});
	clustering.ForPixelClusters([this,&labels](unsigned int cid, const dasp::Cluster& c, unsigned int pid, const dasp::Point& p) {
		labels[pid] = cluster_labels[cid];
	});
	// plot segment boundaries
	boundaries_wt = slimage::Image1ub(clustering.width(), clustering.height(), slimage::Pixel1ub{0});
	dasp::plots::PlotEdges(boundaries_wt, labels, slimage::Pixel1ub{255}, 1);
}

void Segmentation::createLabelsFromBoundaries(const Superpixels& clusters, float threshold)
{
	struct Vertex {
		 // superpixel id
		unsigned int sid;
	};
	struct Edge {
		// edge weight
		float weight;
	};
	// create a graph will all edges with costs >= threshold
	typedef boost::adjacency_list<boost::listS, boost::vecS, boost::undirectedS, Vertex, Edge> MyGraph;
	MyGraph graph;
	std::vector<MyGraph::vertex_descriptor> vids(segmentation_graph.numNodes());
	for(unsigned int i=0; i<segmentation_graph.numNodes(); i++) {
		MyGraph::vertex_descriptor vid = boost::add_vertex(graph);
		graph[vid].sid = i;
		vids[i] = vid;
	}
	for(const graph::Edge& e : segmentation_graph.getEdges()) {
		// take only edges with weight < threshold
		if(e.cost > threshold) {
			continue;
		}
		MyGraph::edge_descriptor eid;
		bool ok;
		boost::tie(eid, ok) = boost::add_edge(vids[e.a], vids[e.b], graph);
		assert(ok);
		graph[eid].weight = e.cost;
	}
	// compute connected components
	typedef std::map<MyGraph::vertex_descriptor, MyGraph::vertices_size_type> component_type;
	component_type component;
	boost::associative_property_map< component_type > component_map(component);
	segment_count = boost::connected_components(graph, component_map);
	cluster_labels.resize(segmentation_graph.numNodes());
	for(unsigned int i=0; i<segmentation_graph.numNodes(); i++) {
		cluster_labels[i] = component[vids[i]];
	}
}

void Segmentation::ucm(const Superpixels& clusters, float threshold)
{
	// every superpixel is one region
	cluster_labels.resize(segmentation_graph.numNodes());
	for(unsigned int i=0; i<cluster_labels.size(); i++) {
		cluster_labels[i] = i;
	}
	// sort edges by weight
	std::vector<graph::Edge> edges = segmentation_graph.getEdges();
	std::sort(edges.begin(), edges.end(), [](const graph::Edge& x, const graph::Edge& y){ return x.cost < y.cost; });
	// cut one edge after another
	for(unsigned int k=0; k<edges.size(); k++) {
		if(edges[k].cost >= threshold) {
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
	relabel();
}

std::vector<slimage::Pixel3ub> Segmentation::computeSegmentColors(const Superpixels& clusters) const
{
//	return plots::CreateRandomColors(segment_count);
	struct Pair { Eigen::Vector3f val; unsigned int cnt; };
	std::vector<Pair> segment_center_sum(segment_count);
	for(unsigned int i=0; i<cluster_labels.size(); i++) {
		Pair& p = segment_center_sum[cluster_labels[i]];
		p.val += clusters.cluster[i].center.color;
		p.cnt ++;
	}
	std::vector<slimage::Pixel3ub> colors(segment_count);
	for(unsigned int i=0; i<colors.size(); i++) {
		Eigen::Vector3f c = segment_center_sum[i].val / static_cast<float>(segment_center_sum[i].cnt);
		c = clusters.ColorToRGB(c);
		slimage::conversion::Convert(slimage::Pixel3f{{c[0],c[1],c[2]}}, colors[i]);
	}
	return colors;
}

slimage::Image3ub Segmentation::computeLabelImage(const Superpixels& clusters) const
{
	std::vector<slimage::Pixel3ub> colors = computeSegmentColors(clusters);
	slimage::Image3ub vis_img(clusters.width(), clusters.height(), slimage::Pixel3ub{{0,0,0}});
	clusters.ForPixelClusters([this,&vis_img,&colors](unsigned int cid, const dasp::Cluster& c, unsigned int pid, const dasp::Point& p) {
		vis_img[pid] = colors[cluster_labels[cid]];
	});
	return vis_img;
}

struct Entry {
	unsigned int a, b;
	float w;
};

typedef float Real;
typedef Eigen::MatrixXf Mat;
typedef Eigen::VectorXf Vec;

void SolveSpectralDense(const std::vector<Entry>& entries, unsigned int n, Vec& ew, Mat& ev)
{
	// creating matrices
	Mat W = Mat::Zero(n,n);
	std::vector<float> Di(n, 0.0f);
	for(const Entry& e : entries) {
		W(e.a, e.b) = e.w;
		W(e.b, e.a) = e.w;
		Di[e.a] += e.w;
		Di[e.b] += e.w;
	}
#ifdef SEGS_VERBOSE
	std::cout << "SpectralSegmentation: W non zero percentage = " << static_cast<float>(2 * entries.size()) / static_cast<float>(n*n) << std::endl;
#endif
	// connect disconnected segments to everything -> ARGH!
	for(unsigned int i=0; i<n; i++) {
		float& di = Di[i];
		if(di == 0) {
#ifdef SEGS_VERBOSE
			std::cout << "Cluster " << i << " has no connections! " << std::endl;
#endif
			// connect the disconnected cluster to all other clusters with a very small weight
			di = 1.0f;
			float q = di / static_cast<float>(n-1);
			for(unsigned int j=0; j<n; j++) {
				if(j == i) continue;
				W(i,j) = q;
				W(j,i) = q;
			}
		}
	}
	// compute matrices D = diagonal(Di) and A = D - W
	Mat D = Mat::Zero(n,n);
	Mat A = -W;
	for(unsigned int i=0; i<n; i++) {
		Real di = Di[i];
		BOOST_ASSERT(di > static_cast<Real>(0));
		A(i,i) += di;
		D(i,i) = di;
	}
	// solve eigensystem
	Eigen::GeneralizedSelfAdjointEigenSolver<Mat> solver;
	solver.compute(A, D); // only need some eigenvectors!
#ifdef SEGS_VERBOSE
	std::cout << "SpectralSegmentation: GeneralizedSelfAdjointEigenSolver says " << solver.info() << std::endl;
#endif
	// return eigenvectors and eigenvalues
	ew = solver.eigenvalues();
	ev = solver.eigenvectors();
}

//void SolveSpectralSparse(const std::vector<Entry>& entries, unsigned int n, Vec& ew, Mat& ev)
//{
//	std::vector<float> Di(n, 0.0f);
//	std::vector<Eigen::Triplet<Real>> W_triplets;
//	W_triplets.reserve(2*entries.size());
//	for(const Entry& e : entries) {
//		Di[e.a] += e.w;
//		Di[e.b] += e.w;
//		W_triplets.push_back(Eigen::Triplet<Real>(e.a,e.b,e.w));
//		W_triplets.push_back(Eigen::Triplet<Real>(e.b,e.a,e.w));
//	}
//	std::vector<Eigen::Triplet<Real>> D_triplets;
//	// connect disconnected segments to everything -> ARGH!
//	for(unsigned int i=0; i<n; i++) {
//		float& di = Di[i];
//		if(di == 0) {
//			std::cout << "Cluster " << i << " has no connections! " << std::endl;
//			// connect the disconnected cluster to all other clusters with a very small weight
//			di = 1.0f;
//			float q = di / static_cast<float>(n-1);
//			for(unsigned int j=0; j<n; j++) {
//				if(j == i) continue;
//				W_triplets.push_back(Eigen::Triplet<Real>(i,j,q));
//				W_triplets.push_back(Eigen::Triplet<Real>(j,i,q));
//			}
//		}
//		D_triplets.push_back(Eigen::Triplet<Real>(i,i,di));
//	}
//	// create matrices
//	Eigen::SparseMatrix<Real> W(n,n);
//	W.setFromTriplets(W_triplets.begin(), W_triplets.end());
//	Eigen::SparseMatrix<Real> D(n,n);
//	W.setFromTriplets(D_triplets.begin(), D_triplets.end());
//	Eigen::SparseMatrix<Real> A = D - W;
//	// solve eigensystem
//	Eigen::GeneralizedSelfAdjointEigenSolver< Eigen::SparseSelfAdjointView<Eigen::SparseMatrix<Real>,1u> > solver;
//	solver.compute(A, D); // only need some eigenvectors!
//	std::cout << "SpectralSegmentation: GeneralizedSelfAdjointEigenSolver says " << solver.info() << std::endl;
//	// return eigenvectors and eigenvalues
//	ew = solver.eigenvalues();
//	ev = solver.eigenvectors();
//}

void SolveSpectral(const std::vector<Entry>& entries, unsigned int n, Vec& ew, Mat& ev)
{
	SolveSpectralDense(entries, n, ew, ev);

	//	{	// DEBUG
	//		std::ofstream ofs_D("/tmp/spectral_D.csv");
	//		std::ofstream ofs_W("/tmp/spectral_W.csv");
	//		for(unsigned int i=0; i<n; i++) {
	//			for(unsigned int j=0; j<n; j++) {
	//				ofs_D << D(i,j);
	//				ofs_W << W(i,j);
	//				if(j+1 == n) {
	//					ofs_D << std::endl;
	//					ofs_W << std::endl;
	//				}
	//				else {
	//					ofs_D << ",";
	//					ofs_W << ",";
	//				}
	//			}
	//		}
	//	}	// DEBUG

}

Segmentation SpectralSegmentation(const Superpixels& clusters, const SpectralSettings& settings)
{
	const bool cOnlyConcaveEdges = true;

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

	const unsigned int cNEV = settings.num_eigenvectors;
	const float cWeightRho = 0.01f; // 640x480 clusters would yield 0.1 which is used in gPb
	unsigned int n = clusters.clusterCount();
#ifdef SEGS_VERBOSE
	std::cout << "SpectralSegmentation: n = " << n << std::endl;
#endif
	// create local neighbourhood graph
	Superpixels::NeighborGraphSettings Gnb_settings;
	Gnb_settings.cut_by_spatial = false;
	Gnb_settings.min_abs_border_overlap = 2;
	Gnb_settings.min_border_overlap = 0.00f;
	Gnb_settings.cost_function = Superpixels::NeighborGraphSettings::SpatialNormalColor;
	graph::Graph Gnb = clusters.CreateNeighborhoodGraph(Gnb_settings);
	BOOST_ASSERT(n == Gnb.nodes_);
	// create W matrix from neighbourhood graph
	Vec edge_connectivity(Gnb.numEdges());
	std::vector<Entry> entries;
	for(unsigned int i=0; i<Gnb.numEdges(); i++) {
		graph::Edge& e = Gnb.getEdges()[i];
		BOOST_ASSERT(e.a != e.b);
		// compute individual edge distances
		float w_maha_color;
		float w_maha_spatial;
		float w_maha_normal;
		// FIXME metric needs central place!
		if(clusters.opt.weight_image == 0.0f) {
			w_maha_color = e.c_color / (std::sqrt(static_cast<float>(clusters.clusterCount())) * cWeightRho);
			w_maha_spatial = std::max(0.0f, std::min(4.0f, e.c_world/4.0f - 1.0f)); // distance of 2 indicates normal distance
			if(cOnlyConcaveEdges) {
				// only use concave edges
				const Point& ca = clusters.cluster[e.a].center;
				const Point& cb = clusters.cluster[e.b].center;
				Eigen::Vector3f d = (cb.world - ca.world).normalized();
				float ca1 = ca.computeNormal().dot(d);
				float ca2 = cb.computeNormal().dot(d);
				float w = ca1 - ca2;
				//float w = std::acos(ca2) - std::acos(ca1);
				if(w < 0.0f) {
					w_maha_normal = 0.0f;
				}
				else {
					w_maha_normal = 3.0f * w;
					//w_maha_normal = 3.0f * (1 - std::cos(w));
				}
			}
			else {
				w_maha_normal = 3.0f*e.c_normal;
			}
		}
		else {
			w_maha_color = 4.0f * e.c_color / (std::sqrt(static_cast<float>(clusters.clusterCount())) * cWeightRho);
			w_maha_spatial = 0.0f;
			w_maha_normal = 0.0f;
		}
		// compute total edge connectivity
//		std::cout << settings.w_spatial << " " << w_maha_spatial << " " << settings.w_color << " " << w_maha_color << " " << settings.w_normal << " " << w_maha_normal << std::endl;
		float w = std::exp(-(settings.w_spatial*w_maha_spatial + settings.w_color*w_maha_color + settings.w_normal*w_maha_normal));
//		float w = std::exp(-w_maha_spatial);
//		float w = std::exp(-w_maha_normal);
//		float w = std::exp(-w_maha_color);
		edge_connectivity[i] = w;
		// write in W and D matrix
		entries.push_back({e.a, e.b, w});
	}

#ifdef SEGS_VERBOSE
	std::cout << "Edve connectivity: min=" << edge_connectivity.minCoeff() << ", max=" << edge_connectivity.maxCoeff() << std::endl;
#endif

	// compute edge border pixels
	std::vector<std::vector<unsigned int>> border_pixels = clusters.ComputeBorderPixels(Gnb);

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
	SolveSpectral(entries, n, result_ew, result_ev);


	unsigned int n_used_ew = std::min(n - 1, cNEV);
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
	Vec edge_weight = Vec::Zero(Gnb.numEdges());
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
		Vec e_k = Vec::Zero(Gnb.numEdges());
		for(unsigned int eid=0; eid<Gnb.numEdges(); eid++) {
			const graph::Edge& e = Gnb.getEdges()[eid];
			e_k[eid] = std::abs(ev[e.a] - ev[e.b]);
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

	// original edge connectivity graph
	graph::Graph graph_original(Gnb.numNodes());
	for(unsigned int eid=0; eid<Gnb.getEdges().size(); eid++) {
		graph::Edge e = Gnb.getEdges()[eid];
		e.cost = edge_connectivity[eid];
		graph_original.add(e);
	}

	// create superpixel neighbourhood graph with edge strength
	graph::Graph graph(Gnb.numNodes());
	for(unsigned int eid=0; eid<Gnb.getEdges().size(); eid++) {
		graph::Edge e = Gnb.getEdges()[eid];
		e.cost = edge_weight[eid];
		graph.add(e);
	}

	Segmentation segs;
	segs.original_graph = graph_original;
	segs.original_boundaries = Segmentation::CreateBorderImageInv(clusters.width(), clusters.height(), graph_original, border_pixels);
	segs.boundaries = Segmentation::CreateBorderImage(clusters.width(), clusters.height(), graph, border_pixels);
	segs.segmentation_graph = graph;

#ifdef SEGS_DBG_SHOWGUI
	slimage::gui::Show("boundaries", segs.boundaries, 0.03f);
#endif

#ifdef SEGS_DBG_SHOWGUI
	slimage::gui::WaitForKeypress();
#endif

	return segs;
}

Segmentation MinCutSegmentation(const Superpixels& clusters)
{
	graph::Graph Gn = clusters.CreateNeighborhoodGraph();
	Segmentation seg;
	// graph segmentation
	seg.segmentation_graph = graph::MinimalSpanningCutting(Gn, clusters.opt.segment_threshold, &seg.cluster_labels);
	// remap labels to get a continuous interval of labels
	seg.relabel();
	return seg;
}

void Segmentation::relabel()
{
	std::set<unsigned int> unique_labels_set(cluster_labels.begin(), cluster_labels.end());
	std::vector<unsigned int> unique_labels(unique_labels_set.begin(), unique_labels_set.end());
	segment_count = unique_labels.size();
	for(unsigned int& x : cluster_labels) {
		auto it = std::find(unique_labels.begin(), unique_labels.end(), x);
		x = it - unique_labels.begin();
	}
}

}
