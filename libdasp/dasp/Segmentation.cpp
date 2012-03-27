/*
 * Segmentation.cpp
 *
 *  Created on: Mar 26, 2012
 *      Author: david
 */

#include "Segmentation.hpp"
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <boost/assert.hpp>
#include <boost/format.hpp>
#include <fstream>
#include <iostream>

namespace dasp
{

std::vector<slimage::Image3ub> cSegmentationDebug;

slimage::Image1f SpectralSegmentation(const Clustering& clusters)
{
	const bool cDebugSaveImages = false;

	const unsigned int cNEV = 48;
	const float cWeightRho = 1.0f;
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
	Eigen::MatrixXf W = Eigen::MatrixXf::Zero(n,n);
	std::vector<float> Di(n, 0.0f);
	for(const graph::Edge& e : Gnb.edges) {
		BOOST_ASSERT(e.a != e.b);
		float w = std::exp(-e.cost / cWeightRho);
		W(e.a, e.b) = w;
		W(e.b, e.a) = w;
		Di[e.a] += w;
		Di[e.b] += w;
	}
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
	Eigen::MatrixXf D = Eigen::MatrixXf::Zero(n,n);
	Eigen::MatrixXf A = -W;
	for(unsigned int i=0; i<n; i++) {
		float di = Di[i];
		BOOST_ASSERT(di > 0.0f);
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
	Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXf> solver;
	solver.compute(A, D); // only need some eigenvectors!
	std::cout << "SpectralSegmentation: GeneralizedSelfAdjointEigenSolver says " << solver.info() << std::endl;
	unsigned int n_used_ew = std::min(n - 1, cNEV);
	std::cout << "Eigenvalues = " << solver.eigenvalues().topRows(n_used_ew + 1).transpose() << std::endl;
	if(cDebugSaveImages) {	// DEBUG
		cSegmentationDebug.clear();
		// create image from eigenvectors (omit first)
		for(unsigned int k=0; k<n_used_ew; k++) {
			// get k-th eigenvector
			Eigen::VectorXf ev = solver.eigenvectors().col(k + 1);
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
		}
	}	// DEBUG
	// compute edge border pixels
	std::vector<std::vector<unsigned int>> border_pixels = clusters.ComputeBorderPixels(Gnb);
	Eigen::VectorXf edge_weight = Eigen::VectorXf::Zero(Gnb.edges.size());
	// later we weight by eigenvalues
	// find a positive eigenvalue (need to do this because of ugly instabilities ...
	float ew_pos = -1.0f;
	for(unsigned int i=0; ; i++) {
		if(solver.eigenvalues()[i] > 0) {
			// FIXME magic to get a not too small eigenvalue
//			unsigned int x = (n_used_ew + i)/2;
			unsigned int x = i + 5;
			ew_pos = solver.eigenvalues()[x];
			break;
		}
	}
	// compute normalized weights from eigenvalues
	Eigen::VectorXf weights = Eigen::VectorXf::Zero(n_used_ew);
	for(unsigned int k=0; k<n_used_ew; k++) {
		float ew = solver.eigenvalues()[k + 1];
		if(ew <= ew_pos) {
			ew = ew_pos;
		}
		weights[k] = 1.0f / std::sqrt(ew);
	}
	std::cout << "Weights = " << weights.transpose() << std::endl;
	// look into first eigenvectors
	for(unsigned int k=0; k<n_used_ew; k++) {
		// weight by eigenvalue
		float w = weights[k];// / w_sum;
		Eigen::VectorXf ev = solver.eigenvectors().col(k + 1);
		// for each edge
		for(unsigned int eid=0; eid<Gnb.edges.size(); eid++) {
			const graph::Edge& e = Gnb.edges[eid];
			// sum difference between eigenvector values for superpixels
			float d = std::abs(ev[e.a] - ev[e.b]);
			edge_weight[eid] += w * d;
		}
	}
	// normalize edge weights to [0,1]
	std::cout << "Edge weights: min=" << edge_weight.minCoeff() << ", max=" << edge_weight.maxCoeff() << std::endl;
	edge_weight /= edge_weight.maxCoeff();
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
			result[pid] = edge_weight[eid];
		}
	}
	return result;
}

slimage::Image1f ComputeBoundary(const Clustering& clusters)
{
	return SpectralSegmentation(clusters);
}

}
