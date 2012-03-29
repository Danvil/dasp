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
#include <Eigen/Eigenvalues>
#include <opencv2/opencv.hpp>
#include <boost/assert.hpp>
#include <boost/format.hpp>
#include <fstream>
#include <iostream>
#include <set>

//#define SEGS_DBG_SHOWGUI
//#define SEGS_DBG_CREATE_EV_IMAGES
//#define SEGS_DBG_PRINT

namespace dasp
{

std::vector<slimage::Image3ub> cSegmentationDebug;

void Segmentation::createBoundariesFromLabels(const Clustering& clustering)
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

Segmentation SpectralSegmentation(const Clustering& clusters)
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

	typedef float Real;
	typedef Eigen::MatrixXf Mat;
	typedef Eigen::VectorXf Vec;

	const unsigned int cNEV = 16;
	const float cWeightRho = 0.01f; // 640x480 clusters would yield 0.1 which is used in gPb
	unsigned int n = clusters.clusterCount();
	std::cout << "SpectralSegmentation: n = " << n << std::endl;
	// create local neighbourhood graph
	Clustering::NeighborGraphSettings Gnb_settings;
	Gnb_settings.cut_by_spatial = false;
	Gnb_settings.min_border_overlap = 0.05f;
	Gnb_settings.cost_function = Clustering::NeighborGraphSettings::SpatialNormalColor;
	graph::Graph Gnb = clusters.CreateNeighborhoodGraph(Gnb_settings);
	BOOST_ASSERT(n == Gnb.nodes_);
	// create W matrix from neighbourhood graph
	Vec edge_connectivity(Gnb.edges.size());
	Mat W = Mat::Zero(n,n);
	std::vector<float> Di(n, 0.0f);
	for(unsigned int i=0; i<Gnb.edges.size(); i++) {
		const graph::Edge& e = Gnb.edges[i];
		BOOST_ASSERT(e.a != e.b);
		// compute individual edge distances
		float w_maha_color = e.c_color / (std::sqrt(static_cast<float>(clusters.clusterCount())) * cWeightRho);
		float w_maha_world = std::max(0.0f, std::min(4.0f, e.c_world/2.0f - 1.0f)); // distance of 2 indicates normal distance
		float w_maha_normal;
		if(cOnlyConcaveEdges) {
			// only use concave edges
			const Point& ca = clusters.cluster[e.a].center;
			const Point& cb = clusters.cluster[e.b].center;
			Eigen::Vector3f d = (cb.world - ca.world).normalized();
			float ca1 = ca.computeNormal().dot(d);
			float ca2 = cb.computeNormal().dot(d);
			float w = ca1 - ca2;
			if(w < 0.0f) {
				w_maha_normal = 0.0f;
			}
			else {
				w_maha_normal = 3.0f*w;
			}
		}
		else {
			w_maha_normal = 3.0f*e.c_normal;
		}
		// compute total edge connectivity
		float w = std::exp(-(w_maha_color + w_maha_world + w_maha_normal));
//		float w = std::exp(-w_maha_world);
//		float w = std::exp(-w_maha_normal);
//		float w = std::exp(-w_maha_color);
		edge_connectivity[i] = w;
		// write in W and D matrix
		W(e.a, e.b) = w;
		W(e.b, e.a) = w;
		Di[e.a] += w;
		Di[e.b] += w;
	}

	std::cout << "Edve connectivity: min=" << edge_connectivity.minCoeff() << ", max=" << edge_connectivity.maxCoeff() << std::endl;

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

	std::cout << "SpectralSegmentation: W non zero percentage = " << static_cast<float>(2 * Gnb.edges.size()) / static_cast<float>(n*n) << std::endl;
	// connect disconnected segments to everything -> ARGH!
	for(unsigned int i=0; i<n; i++) {
		float& di = Di[i];
		if(di == 0) {
			std::cout << "Cluster " << i << " has no connections! pixels=" << clusters.cluster[i].pixel_ids.size() << std::endl;
			di = 1.0f; // FIXME wtf!
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
	// solve eigensystem
	Eigen::GeneralizedSelfAdjointEigenSolver<Mat> solver;
	solver.compute(A, D); // only need some eigenvectors!
	std::cout << "SpectralSegmentation: GeneralizedSelfAdjointEigenSolver says " << solver.info() << std::endl;
	unsigned int n_used_ew = std::min(n - 1, cNEV);
	std::cout << "Eigenvalues = " << solver.eigenvalues().topRows(n_used_ew + 1).transpose() << std::endl;
#ifdef SEGS_DBG_CREATE_EV_IMAGES
	{
		cSegmentationDebug.clear();
		// create image from eigenvectors (omit first)
		for(unsigned int k=0; k<std::min(cNEV,3u); k++) {
			// get k-th eigenvector
			Vec ev = solver.eigenvectors().col(k + 1);
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
	Vec edge_weight = Vec::Zero(Gnb.edges.size());
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
		Real ew = solver.eigenvalues()[k + 1];
		if(ew <= Real(0)) {
			// omit if eigenvalue is not positive
			continue;
		}
		float w = 1.0f / std::sqrt(ew);
		// get eigenvector and normalize
		Vec ev = solver.eigenvectors().col(k + 1);
		ev = (ev - ev.minCoeff()*Vec::Ones(ev.rows())) / (ev.maxCoeff() - ev.minCoeff());
		// for each edge compute difference of eigenvector values
		Vec e_k = Vec::Zero(Gnb.edges.size());
		for(unsigned int eid=0; eid<Gnb.edges.size(); eid++) {
			const graph::Edge& e = Gnb.edges[eid];
			e_k[eid] = std::abs(ev[e.a] - ev[e.b]);
		}
		std::cout << "w=" << w << " e_k.maxCoeff()=" << e_k.maxCoeff() << std::endl;
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


	// normalize edge weights to [0,1]
	std::cout << "Edge weights: min=" << edge_weight.minCoeff() << ", max=" << edge_weight.maxCoeff() << std::endl;
	std::cout << "Edge weights: median=" << edge_weight[edge_weight.rows()/2] << std::endl;

//	edge_weight /= edge_weight[(95*edge_weight.rows())/100];
//	edge_weight /= edge_weight.maxCoeff();
//	std::cout << "Edge weights = " << edge_weight.transpose() << std::endl;

	// write value to border pixels
	slimage::Image1f result(clusters.width(), clusters.height(), slimage::Pixel1f{0.0f});
//	for(unsigned int y=1; y<clusters.height()-1; y++) {
//		for(unsigned int x=1; x<clusters.width()-1; x++) {
//			if(clusters.points(x,y).isInvalid()) {
//				continue;
//			}
//			if(clusters.points(x-1,y).isInvalid()
//			|| clusters.points(x+1,y).isInvalid()
//			|| clusters.points(x,y-1).isInvalid()
//			|| clusters.points(x,y+1).isInvalid()) {
//				result(x,y) = 100.0f;
//			}
//		}
//	}
	for(unsigned int eid=0; eid<Gnb.edges.size(); eid++) {
		for(unsigned int pid : border_pixels[eid]) {
			result[pid] = static_cast<float>(edge_weight[eid]);
		}
	}

#ifdef SEGS_DBG_SHOWGUI
	slimage::gui::Show("segs", result, 0.01f);
	slimage::gui::WaitForKeypress();
#endif

	// FIXME create and segmentation graph labels!

	Segmentation segs;
	segs.boundaries = result;
	return segs;
}

Segmentation MinCutSegmentation(const Clustering& clusters)
{
	graph::Graph Gn = clusters.CreateNeighborhoodGraph();
	Segmentation seg;
	// graph segmentation
	seg.segmentation_graph = graph::MinimalSpanningCutting(Gn, clusters.opt.segment_threshold, &seg.cluster_labels);
	// remap labels to get a continuous interval of labels
	std::set<unsigned int> unique_labels_set(seg.cluster_labels.begin(), seg.cluster_labels.end());
	std::vector<unsigned int> unique_labels(unique_labels_set.begin(), unique_labels_set.end());
	seg.segment_count = unique_labels.size();
	for(unsigned int& x : seg.cluster_labels) {
		auto it = std::find(unique_labels.begin(), unique_labels.end(), x);
		x = it - unique_labels.begin();
	}
	return seg;
}



}
