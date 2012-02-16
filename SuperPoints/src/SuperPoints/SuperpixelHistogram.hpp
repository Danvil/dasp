/*
 * SuperpixelHistogram.hpp
 *
 *  Created on: Feb 15, 2012
 *      Author: david
 */

#ifndef SUPERPIXELHISTOGRAM_HPP_
#define SUPERPIXELHISTOGRAM_HPP_
//----------------------------------------------------------------------------//
#include <Danvil/Statistics/Histogram.hpp>
#include <eigen3/Eigen/Dense>
#include <vector>
#include <cmath>
//----------------------------------------------------------------------------//
namespace dasp {
//----------------------------------------------------------------------------//

struct SuperpixelState
{
	Eigen::Vector3f color;
	Eigen::Vector3f normal;
};

struct SuperpixelHistogram
{
private:
	/** color value */
	SuperpixelState state_;

	/** feature vector
	 * difference in color chroma
	 * difference in color intensity
	 * angle between normals
	 */
public:
	Eigen::VectorXf hist_chroma_, hist_intensity_, hist_normal_;

	static const unsigned int cBucketCountChroma = 32;
	static const unsigned int cBucketCountIntensity = 24;
	static const unsigned int cBucketCountNormal = 16;

public:
	SuperpixelHistogram()
	{
		state_.color = Eigen::Vector3f::Zero();
		state_.normal = Eigen::Vector3f::UnitZ();
		clear();
	}

	SuperpixelHistogram(const SuperpixelState& state)
	: state_(state) {
		clear();
	}

	void clear() {

		hist_chroma_ = Eigen::VectorXf::Zero(cBucketCountChroma);
		hist_intensity_ = Eigen::VectorXf::Zero(cBucketCountIntensity);
		hist_normal_ = Eigen::VectorXf::Zero(cBucketCountNormal);
	}

	void add(const SuperpixelState& x) {
		Eigen::Vector3f f = ComputeFeature(x);
		addFeature(f);
	}

	static float Distance(const SuperpixelHistogram& x, const SuperpixelHistogram& y) {
		static Eigen::MatrixXf AC = ComputeSimilarityChroma();
		static Eigen::MatrixXf AI = ComputeSimilarityIntensity();
		static Eigen::MatrixXf AN = ComputeSimilarityNormal();
		// rgb distance
		float d_rgb = (x.state_.color - y.state_.color).norm();
		float d_C = FeatureDistance(x.hist_chroma_, y.hist_chroma_, AC);
		float d_I = FeatureDistance(x.hist_intensity_, y.hist_intensity_, AI);
		float d_N = FeatureDistance(x.hist_normal_, y.hist_normal_, AN);
		return d_rgb + 0.1f * (d_C + d_I + d_N);
	}

private:
	Eigen::Vector3f ComputeFeature(const SuperpixelState& x) const {
		// compute delta chroma (angle between colors)
		float dC = state_.color.dot(x.color) / (state_.color.norm() * x.color.norm());
		// compute delta intensity (difference between color intensities)
		float dI = std::abs(state_.color.sum() - x.color.sum());
		// compute delta normal (angle between normals)
		float dN = state_.normal.dot(x.normal); // normals have unit length
		Eigen::Vector3f f;
		f << dC, dI, dN;
		return f;
	}

	void addFeature(const Eigen::Vector3f& f) {
		float fC = f[0];
		float fI = f[1];
		float fN = f[2];
		addChroma(fC);
		addIntensity(fI);
		addNormal(fN);
	}

	void addChroma(float v) {
		unsigned int b = FindBucket(v, 0.0f, 1.0f, cBucketCountChroma);
		hist_chroma_[b] += 1.0f;
	}

	void addIntensity(float v) {
		unsigned int b = FindBucket(v, 0.0f, 1.0f, cBucketCountIntensity);
		hist_intensity_[b] += 1.0f;
	}

	void addNormal(float v) {
		const float cPi = 3.141592654f;
		unsigned int b = FindBucket(v, -cPi, +cPi, cBucketCountNormal);
		hist_normal_[b] += 1.0f;
	}

	static unsigned int FindBucket(float v, float min, float max, unsigned int buckets) {
		int b = static_cast<unsigned int>(std::floor(float(v - min) * float(buckets) / float(max - min)));
		if(b < 0) {
			return 0;
		}
		else if(b >= int(buckets)) {
			return buckets - 1;
		}
		else {
			return b;
		}
	}

	static Eigen::MatrixXf ComputeSimilarityChroma() {
		Eigen::MatrixXf A = Eigen::MatrixXf::Zero(cBucketCountChroma, cBucketCountChroma);
		for(int i=0; i<int(cBucketCountChroma); i++) {
			for(int j=0; j<int(cBucketCountChroma); j++) {
				// use cluster distance
				float d = static_cast<float>(std::abs(i - j)) / static_cast<float>(cBucketCountChroma);
				A(i,j) = std::max(0.0f, 1.0f - d);
			}
		}
		return A;
	}

	static Eigen::MatrixXf ComputeSimilarityIntensity() {
		Eigen::MatrixXf A = Eigen::MatrixXf::Zero(cBucketCountIntensity, cBucketCountIntensity);
		for(int i=0; i<int(cBucketCountIntensity); i++) {
			for(int j=0; j<int(cBucketCountIntensity); j++) {
				// use cluster distance
				float d = static_cast<float>(std::abs(i - j)) / static_cast<float>(cBucketCountIntensity);
				A(i,j) = std::max(0.0f, 1.0f - d);
			}
		}
		return A;
	}

	static Eigen::MatrixXf ComputeSimilarityNormal() {
		Eigen::MatrixXf A = Eigen::MatrixXf::Zero(cBucketCountNormal, cBucketCountNormal);
		for(int i=0; i<int(cBucketCountNormal); i++) {
			for(int j=0; j<int(cBucketCountNormal); j++) {
				// use cluster distance
				// use cluster distance
				int u = std::abs(i - j);
				// but 0 and 2pi are the same!
				int uu = int(cBucketCountNormal) - u;
				if(uu < u) {
					u = uu;
				}
				// max distance is half!
				float d = static_cast<float>(u) / static_cast<float>(cBucketCountNormal) * 2.0f;
				A(i,j) = std::max(0.0f, 1.0f - d);
			}
		}
		return A;
	}

//	static Eigen::MatrixXf DensityToSimilarity(const Eigen::MatrixXf& D) {
//		Eigen::MatrixXf A = D;
//		float m = A.maxCoeff();
//		BOOST_ASSERT(m > 0.0f);
//		return Eigen::MatrixXf::Ones() - A / m;
//	}

	static float FeatureDistance(const Eigen::VectorXf& x, const Eigen::VectorXf& y, const Eigen::MatrixXf& A) {
		/// For more details on this distance see the paper:
		///  The Quadratic-Chi Histogram Distance Family
		///  Ofir Pele, Michael Werman
		///  ECCV 2010
		const float m = 0.5f;
		Eigen::VectorXf z = A * (x + y);
		Eigen::VectorXf d = x - y;
		for(int i=0; i<z.rows(); i++) {
			float u = z[i];
			if(u == 0.0f) {
				u = 1.0f;
			}
			else {
				u = std::pow(u, m);
			}
			d[i] /= u;
		}
		float d2 = d.transpose() * A * d;
		return std::sqrt(std::max(0.0f, d2));
	}

};

//----------------------------------------------------------------------------//
}
//----------------------------------------------------------------------------//
#endif
