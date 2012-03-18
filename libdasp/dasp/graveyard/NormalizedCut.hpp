/*
 * NormalizedCut.hpp
 *
 *  Created on: Feb 2, 2012
 *      Author: david
 */

#ifndef NORMALIZEDCUT_HPP_
#define NORMALIZEDCUT_HPP_

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Eigenvalues>
#include <vector>
#include <algorithm>

namespace NormalizedCut
{
	std::vector<unsigned int> NCut(const Eigen::MatrixXf& W) {
		assert(W.rows() == W.cols());
		unsigned int n = W.rows();
		Eigen::MatrixXf D = Eigen::MatrixXf::Identity(n,n);
		for(unsigned int i=0; i<n; i++) {
			D(i,i) = W.col(i).sum(); // or row?
		}
		Eigen::MatrixXf A = D - W;
		Eigen::SelfAdjointEigenSolver<Eigen::MatrixXf> esolver;
		esolver.compute(A);
		// Find the eigenvector with the second smallest eigenvalue
		Eigen::VectorXf ew = esolver.eigenvalues();
		std::vector<std::pair<unsigned int,float>> ew_vals(n);
		for(unsigned int i=0; i<n; i++) {
			ew_vals[i].first = i;
			ew_vals[i].second = ew[i];
		}
		std::sort(ew_vals.begin(), ew_vals.end(), [](const std::pair<unsigned int,float>& x, const std::pair<unsigned int,float>& y) { return x.second < y.second; });
		unsigned int ew_min = ew_vals[1].first; // second smallest!!
		// Use the eigenvector with the second smallest eigenvalue to bipartition the graph (e.g. grouping according to sign).
		ew.minCoeff(&ew_min);
		Eigen::VectorXf v = esolver.eigenvectors().col(ew_min);
		std::vector<float> v_vals(n);
		for(unsigned int i=0; i<n; i++) {
			v_vals[i] = v[i];
		}
		std::sort(v_vals.begin(), v_vals.end());
		float splitpoint = v_vals[v_vals.size() / 2];
		std::vector<unsigned int> groupA;
		groupA.reserve(n);
		for(unsigned int i=0; i<n; i++) {
			if(v[i] > splitpoint) {
				groupA.push_back(i);
			}
		}
		return groupA;
	}
}

#endif /* NORMALIZEDCUT_HPP_ */
