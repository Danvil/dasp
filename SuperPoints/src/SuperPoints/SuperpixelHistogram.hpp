/*
 * SuperpixelHistogram.hpp
 *
 *  Created on: Feb 15, 2012
 *      Author: david
 */

#ifndef SUPERPIXELHISTOGRAM_HPP_
#define SUPERPIXELHISTOGRAM_HPP_
//----------------------------------------------------------------------------//
#include "Histogram.hpp"
#include <Danvil/LinAlg/Eigen.hpp>
#include <Danvil/Statistics/KMeans.hpp>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Eigenvalues>
#include <vector>
#include <cmath>
//----------------------------------------------------------------------------//
namespace dasp {
//----------------------------------------------------------------------------//

struct SuperpixelState
{
	unsigned int x, y;
	Eigen::Vector3f position;
	Eigen::Vector3f color;
	Eigen::Vector3f normal;
	float scala;
};

struct SuperpixelNeighbourhoodGraph
{
	SuperpixelState center_;
	std::vector<SuperpixelState> neighbours_;
};

struct SuperpixelGraph
{
	std::vector<SuperpixelState> nodes_;

	std::vector<std::vector<std::size_t>> node_connections_;

	std::size_t size() const {
		return nodes_.size();
	}

	void createConnections(float threshold) {
		std::size_t n = size();
		node_connections_.resize(n);
		for(std::size_t i=0; i<n; i++) {
			for(std::size_t j=i+1; j<n; j++) {
				float d = (nodes_[i].position - nodes_[j].position).norm();
				// only connect if distance is smaller than threshold
				if(d < threshold) {
					node_connections_[i].push_back(j);
					node_connections_[j].push_back(i);
				}
			}
		}
	}

	SuperpixelNeighbourhoodGraph createNeighbourhoodGraph(unsigned int i) const {
		SuperpixelNeighbourhoodGraph ng;
		ng.center_ = nodes_[i];
		for(std::size_t j : node_connections_[i]) {
			ng.neighbours_.push_back(nodes_[j]);
		}
		return ng;
	}

};

//----------------------------------------------------------------------------//

class ISuperpixelModel
{
public:
	virtual ~ISuperpixelModel() {}

	virtual void train(const SuperpixelGraph& graph, const std::vector<unsigned int>& selection) = 0;

	virtual float evaluateNeighbourhood(const SuperpixelNeighbourhoodGraph& x) const = 0;

	virtual unsigned int labelNeighbourhood(const SuperpixelNeighbourhoodGraph& x) const = 0;

	std::vector<float> evaluate(const SuperpixelGraph& g) const {
		std::vector<float> result(g.size());
		for(std::size_t i=0; i<g.size(); i++) {
			SuperpixelNeighbourhoodGraph ng = g.createNeighbourhoodGraph(i);
			result[i] = evaluateNeighbourhood(ng);
			if(i % 37 == 0) {
				std::cout << i << ": n=" << ng.neighbours_.size() << ", v=" << result[i] << std::endl;
			}
		}
		return result;
	}

	std::vector<unsigned int> label(const SuperpixelGraph& g) const {
		std::vector<unsigned int> result(g.size());
		for(std::size_t i=0; i<g.size(); i++) {
			SuperpixelNeighbourhoodGraph ng = g.createNeighbourhoodGraph(i);
			result[i] = labelNeighbourhood(ng);
		}
		return result;
	}

};

//----------------------------------------------------------------------------//

//struct SuperpixelHistogram
//{
//private:
//	/** color value */
//	SuperpixelState state_;
//
//	/** feature vector
//	 * difference in color chroma
//	 * difference in color intensity
//	 * angle between normals
//	 */
//public:
//	Eigen::VectorXf hist_chroma_, hist_intensity_, hist_normal_;
//
//	static const unsigned int cBucketCountChroma = 32;
//	static const unsigned int cBucketCountIntensity = 24;
//	static const unsigned int cBucketCountNormal = 16;
//
//public:
//	SuperpixelHistogram()
//	{
//		state_.color = Eigen::Vector3f::Zero();
//		state_.normal = Eigen::Vector3f::UnitZ();
//		clear();
//	}
//
//	SuperpixelHistogram(const SuperpixelState& state)
//	: state_(state) {
//		clear();
//	}
//
//	void clear() {
//
//		hist_chroma_ = Eigen::VectorXf::Zero(cBucketCountChroma);
//		hist_intensity_ = Eigen::VectorXf::Zero(cBucketCountIntensity);
//		hist_normal_ = Eigen::VectorXf::Zero(cBucketCountNormal);
//	}
//
//	void add(const SuperpixelState& x) {
//		Eigen::Vector3f f = ComputeFeature(x);
//		addFeature(f);
//	}
//
//	static float Distance(const SuperpixelHistogram& x, const SuperpixelHistogram& y) {
//		static Eigen::MatrixXf AC = ComputeSimilarityChroma();
//		static Eigen::MatrixXf AI = ComputeSimilarityIntensity();
//		static Eigen::MatrixXf AN = ComputeSimilarityNormal();
//		// rgb distance
//		float d_rgb = (x.state_.color - y.state_.color).norm();
//		float d_C = HistogramDistance(x.hist_chroma_, y.hist_chroma_, AC);
//		float d_I = HistogramDistance(x.hist_intensity_, y.hist_intensity_, AI);
//		float d_N = HistogramDistance(x.hist_normal_, y.hist_normal_, AN);
//		return d_rgb + 0.1f * (d_C + d_I + d_N);
//	}
//
//private:
//	Eigen::Vector3f ComputeFeature(const SuperpixelState& x) const {
//		// compute delta chroma (angle between colors)
//		float dC = state_.color.dot(x.color) / (state_.color.norm() * x.color.norm());
//		// compute delta intensity (difference between color intensities)
//		float dI = std::abs(state_.color.sum() - x.color.sum());
//		// compute delta normal (angle between normals)
//		float dN = state_.normal.dot(x.normal); // normals have unit length
//		Eigen::Vector3f f;
//		f << dC, dI, dN;
//		return f;
//	}
//
//	void addFeature(const Eigen::Vector3f& f) {
//		float fC = f[0];
//		float fI = f[1];
//		float fN = f[2];
//		addChroma(fC);
//		addIntensity(fI);
//		addNormal(fN);
//	}
//
//	void addChroma(float v) {
//		unsigned int b = FindBucket(v, 0.0f, 1.0f, cBucketCountChroma);
//		hist_chroma_[b] += 1.0f;
//	}
//
//	void addIntensity(float v) {
//		unsigned int b = FindBucket(v, 0.0f, 1.0f, cBucketCountIntensity);
//		hist_intensity_[b] += 1.0f;
//	}
//
//	void addNormal(float v) {
//		const float cPi = 3.141592654f;
//		unsigned int b = FindBucket(v, -cPi, +cPi, cBucketCountNormal);
//		hist_normal_[b] += 1.0f;
//	}
//
//	static unsigned int FindBucket(float v, float min, float max, unsigned int buckets) {
//		int b = static_cast<unsigned int>(std::floor(float(v - min) * float(buckets) / float(max - min)));
//		if(b < 0) {
//			return 0;
//		}
//		else if(b >= int(buckets)) {
//			return buckets - 1;
//		}
//		else {
//			return b;
//		}
//	}
//
//	static Eigen::MatrixXf ComputeSimilarityChroma() {
//		Eigen::MatrixXf A = Eigen::MatrixXf::Zero(cBucketCountChroma, cBucketCountChroma);
//		for(int i=0; i<int(cBucketCountChroma); i++) {
//			for(int j=0; j<int(cBucketCountChroma); j++) {
//				// use cluster distance
//				float d = static_cast<float>(std::abs(i - j)) / static_cast<float>(cBucketCountChroma);
//				A(i,j) = std::max(0.0f, 1.0f - d);
//			}
//		}
//		return A;
//	}
//
//	static Eigen::MatrixXf ComputeSimilarityIntensity() {
//		Eigen::MatrixXf A = Eigen::MatrixXf::Zero(cBucketCountIntensity, cBucketCountIntensity);
//		for(int i=0; i<int(cBucketCountIntensity); i++) {
//			for(int j=0; j<int(cBucketCountIntensity); j++) {
//				// use cluster distance
//				float d = static_cast<float>(std::abs(i - j)) / static_cast<float>(cBucketCountIntensity);
//				A(i,j) = std::max(0.0f, 1.0f - d);
//			}
//		}
//		return A;
//	}
//
//	static Eigen::MatrixXf ComputeSimilarityNormal() {
//		Eigen::MatrixXf A = Eigen::MatrixXf::Zero(cBucketCountNormal, cBucketCountNormal);
//		for(int i=0; i<int(cBucketCountNormal); i++) {
//			for(int j=0; j<int(cBucketCountNormal); j++) {
//				// use cluster distance
//				// use cluster distance
//				int u = std::abs(i - j);
//				// but 0 and 2pi are the same!
//				int uu = int(cBucketCountNormal) - u;
//				if(uu < u) {
//					u = uu;
//				}
//				// max distance is half!
//				float d = static_cast<float>(u) / static_cast<float>(cBucketCountNormal) * 2.0f;
//				A(i,j) = std::max(0.0f, 1.0f - d);
//			}
//		}
//		return A;
//	}
//
////	static Eigen::MatrixXf DensityToSimilarity(const Eigen::MatrixXf& D) {
////		Eigen::MatrixXf A = D;
////		float m = A.maxCoeff();
////		BOOST_ASSERT(m > 0.0f);
////		return Eigen::MatrixXf::Ones() - A / m;
////	}
//
//};

//----------------------------------------------------------------------------//

class SuperpixelHistogramModel
: public ISuperpixelModel
{
public:
	void train(const SuperpixelGraph& graph, const std::vector<unsigned int>& selection)
	{
		const unsigned int cBins = 15;
		assert(graph.nodes_.size() >= cBins);
		// train model
		std::vector<Danvil::ctLinAlg::Vec3ub> gmm_training;
		gmm_training.reserve(graph.size());
		for(std::size_t i=0; i<graph.size(); i++) {
			 Eigen::VectorXf c = 255.0f * graph.nodes_[i].color;
			 gmm_training.push_back(Danvil::ctLinAlg::Vec3ub(c[0], c[1], c[2]));
		}
		clusters_ = Danvil::Clustering::KMeans(gmm_training, cBins);

		// assure minimal sigma for cluster gaussians
		const float cMinimalSigma = 0.05f;
		for(unsigned int i=0; i<clusters_.cluster_.size(); i++) {
			Eigen::Matrix3f s = Danvil::ctLinAlg::Convert(clusters_.cluster_[i].gaussian().sigma());
			Eigen::SelfAdjointEigenSolver<Eigen::Matrix3f> solver(s);
			Eigen::Vector3f ew = solver.eigenvalues();
			for(unsigned int k=0; k<3; k++) {
				ew[k] = std::max(ew[k], cMinimalSigma*cMinimalSigma);
			}
			s = solver.eigenvectors().inverse() * Eigen::DiagonalMatrix<float,3,3>(ew) * solver.eigenvectors();
			clusters_.cluster_[i].gaussian().SetCovariance(Danvil::ctLinAlg::Convert(s));
		}

		// compute histogram bin similarity
		bin_similarity_ = Eigen::MatrixXf::Zero(cBins, cBins);
		for(unsigned int i=0; i<cBins; i++) {
			bin_similarity_(i,i) = 1.0f;
			for(unsigned int j=i+1; j<cBins; j++) {
				float d = (clusters_.cluster_[i].gaussian().mean() - clusters_.cluster_[j].gaussian().mean()).length();
				bin_similarity_(i,j) = bin_similarity_(j,i) = DistanceToSimilarity(d);
			}
		}

		// DEBUG
		for(unsigned int i=0; i<cBins; i++) {
			auto u = clusters_.cluster_[i].gaussian().mean();
			Eigen::Vector3f v;
			v << u[0], u[1], u[2];
			std::cout << "Gaussian " << i << ": " << v.transpose() << std::endl;
		}
		std::cout << "Bin Similarity" << std::endl;
		std::cout << bin_similarity_ << std::endl;

		// create histograms
		model_hist_.clear();
		for(std::size_t i : selection) {
			SuperpixelNeighbourhoodGraph ng = graph.createNeighbourhoodGraph(i);
			Eigen::VectorXf h = createHistogram(ng);
			std::cout << "hist " << i << ": n=" << ng.neighbours_.size() << ", h=" << h.transpose() << std::endl;
			model_hist_.push_back(h);
		}
	}

	float evaluateNeighbourhood(const SuperpixelNeighbourhoodGraph& x) const {
		Eigen::VectorXf xh = createHistogram(x);
		float c_min = 1e9;
		for(const Eigen::VectorXf& mh : model_hist_) {
			float c = HistogramDistance(xh, mh, bin_similarity_);
			c_min = std::min(c, c_min);
		}
//		return std::max(0.0f, 1.0f - c_min);
		if (c_min < 3.14)
			return 0.5 * cos(c_min) + 0.5;

		return 0.0f;
	}

	Eigen::VectorXf createHistogram(const SuperpixelNeighbourhoodGraph& g) const {
		Eigen::VectorXf h = createHistogram(g.center_.color);
		for(std::size_t i=0; i<g.neighbours_.size(); i++) {
			float d = 1.0f / (1.0f + (g.neighbours_[i].position - g.center_.position).squaredNorm());
			Eigen::VectorXf t = createHistogram(g.neighbours_[i].color);
			h += d * t;
		}
		return h;
	}

	unsigned int labelNeighbourhood(const SuperpixelNeighbourhoodGraph& x) const {
		Danvil::ctLinAlg::Vec3ub cub(255.0f*x.center_.color[0], 255.0f*x.center_.color[1], 255.0f*x.center_.color[2]);
		return clusters_.NearestCluster(cub);
	}

	std::vector<Eigen::Vector3f> getClusterColors() const {
		std::vector<Eigen::Vector3f> colors(clusters_.cluster_.size());
		for(unsigned int i=0; i<clusters_.cluster_.size(); i++) {
			colors[i] = Danvil::ctLinAlg::Convert(clusters_.cluster_[i].gaussian().mean());
		}
		return colors;
	}

private:
	Eigen::VectorXf createHistogram(const Eigen::Vector3f& color) const {
		Danvil::ctLinAlg::Vec3ub cub(255.0f*color[0], 255.0f*color[1], 255.0f*color[2]);
//		std::vector<float> d = clusters_.ClusterProbabilities(Danvil::ctLinAlg::Vec3ub(c[0], c[1], c[2]));
		Eigen::VectorXf h = Eigen::VectorXf::Zero(clusters_.cluster_.size());
//		for(unsigned int i=0; i<clusters_.cluster_.size(); i++) {
//			float p = clusters_.cluster_[i].gaussian().ProbabilityUnscaled(cub);
//			if(p > 0.75f) { // FIXME constant!!!
//				h[i] = 1.0f;
//			}
//		}
		int best = clusters_.NearestCluster(cub);
		h[best] = 1.0f;
		return h;
	}

	inline
	static float DistanceToSimilarity(float d) {
		float u = d / 0.05f;
		return std::exp(-0.5f*u*u);
	}

private:
	Danvil::Clustering::ClusterGroup clusters_;

	Eigen::MatrixXf bin_similarity_;

	std::vector<Eigen::VectorXf> model_hist_;

};

//----------------------------------------------------------------------------//
}
//----------------------------------------------------------------------------//
#endif
