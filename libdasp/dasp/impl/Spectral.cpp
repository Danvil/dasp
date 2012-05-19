/*
 * Spectral.cpp
 *
 *  Created on: May 19, 2012
 *      Author: david
 */

#define SEGS_VERBOSE

#include "Spectral.hpp"
#include <Eigen/Eigenvalues>
#include <boost/assert.hpp>
#ifdef SEGS_VERBOSE
	#include <iostream>
#endif

namespace dasp
{
namespace detail
{

	typedef float Real;
	typedef Eigen::MatrixXf Mat;
	typedef Eigen::VectorXf Vec;

	struct EigenComponent {
		float eigenvalue;
		Vec eigenvector;
	};

	typedef std::vector<EigenComponent> PartialEigenSolution;

	PartialEigenSolution SolveDense(const SpectralGraph& graph, unsigned int num_ev)
	{
		unsigned int dim = boost::num_vertices(graph);
#ifdef SEGS_VERBOSE
		unsigned int num_edges = boost::num_edges(graph);
		std::cout << "SpectralSegmentation: dimension=" << dim
				<< ", num_edges=" << num_edges
				<< ", matrix non-zero elements = " << 100*static_cast<float>(2 * num_edges) / static_cast<float>(dim*dim) << "%" << std::endl;
#endif
		// creating matrices
		Mat W = Mat::Zero(dim,dim);
		std::vector<float> Di(dim, 0.0f);
		for(auto eid : as_range(boost::edges(graph))) {
			unsigned int ea = boost::source(eid, graph);
			unsigned int eb = boost::target(eid, graph);
			float ew = boost::get(boost::edge_weight, graph, eid);
			W(ea, eb) = ew;
			W(eb, ea) = ew;
			Di[ea] += ew;
			Di[eb] += ew;
		}
		// connect disconnected segments to everything
		// FIXME why is this necessary?
		for(unsigned int i=0; i<dim; i++) {
			float& di = Di[i];
			if(di == 0) {
#ifdef SEGS_VERBOSE
				std::cout << "Node " << i << " has no connections! " << std::endl;
#endif
				// connect the disconnected cluster to all other clusters with a very small weight
				di = 1.0f;
				float q = di / static_cast<float>(dim-1);
				for(unsigned int j=0; j<dim; j++) {
					if(j == i) continue;
					W(i,j) = q;
					W(j,i) = q;
				}
			}
		}
		// compute matrices D = diagonal(Di) and A = D - W
		Mat D = Mat::Zero(dim,dim);
		Mat A = -W;
		for(unsigned int i=0; i<dim; i++) {
			Real di = Di[i];
			BOOST_ASSERT(di > static_cast<Real>(0));
			A(i,i) += di;
			D(i,i) = di;
		}
		// solve eigensystem
		Eigen::GeneralizedSelfAdjointEigenSolver<Mat> solver;
		solver.compute(A, D);
#ifdef SEGS_VERBOSE
		std::cout << "SpectralSegmentation: GeneralizedSelfAdjointEigenSolver says " << solver.info() << std::endl;
#endif
		// return eigenvectors and eigenvalues
		PartialEigenSolution solution(std::min(dim, num_ev));
		for(unsigned int i=0; i<solution.size(); i++) {
			solution[i].eigenvalue = solver.eigenvalues()[i];
#ifdef SEGS_VERBOSE
			std::cout << "SpectralSegmentation: eigenvalue #" << i << "=" << solution[i].eigenvalue << std::endl;
#endif
			solution[i].eigenvector = solver.eigenvectors().col(i);
		}
		return solution;
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

	/** Solves the spectral graph theory eigenvalue problem and smallest eigenvalues */
	PartialEigenSolution SpectralPartialEigenSolve(const SpectralGraph& graph, unsigned int num_ev)
	{
		return SolveDense(graph, num_ev);

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

	/** Assembles edge weights from eigenvalues and eigenvectors */
	Vec AssembleEdgeWeights(const SpectralGraph& graph, const PartialEigenSolution& solution)
	{
		Vec edge_weight = Vec::Zero(boost::num_edges(graph));
	//	// later we weight by eigenvalues
	//	// find a positive eigenvalue (need to do this because of ugly instabilities ...
	//	Real ew_pos = -1.0f;
	//	for(unsigned int i=0; ; i++) {
	//		if(solver.eigenvalues()[i] > 0) {
	//			// F IXME magic to get a not too small eigenvalue
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
		// skip first component
		for(unsigned int k=1; k<solution.size(); k++) {
			const EigenComponent& eigen = solution[k];
			// omit if eigenvalue is not positive
			Real ew = eigen.eigenvalue;
			// FIXME this is due to numerical instabilities
			if(ew <= Real(0)) {
				continue;
			}
			// weight by eigenvalue
			float w = 1.0f / std::sqrt(ew);
			// get eigenvector and normalize
			Vec ev = eigen.eigenvector;
			ev = (ev - ev.minCoeff() * Vec::Ones(ev.rows())) / (ev.maxCoeff() - ev.minCoeff());
			// for each edge compute difference of eigenvector values
			Vec e_k = Vec::Zero(edge_weight.rows());
			// FIXME proper edge indexing
			unsigned int eid_index = 0;
			for(auto eid : as_range(boost::edges(graph))) {
				e_k[eid_index] = std::abs(ev[boost::source(eid, graph)] - ev[boost::target(eid, graph)]);
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
		return edge_weight;
	}

}

SpectralGraph SolveSpectral(const SpectralGraph& graph, unsigned int num_ev)
{
	// pick one more because the first one is omitted
	detail::PartialEigenSolution solution = detail::SpectralPartialEigenSolve(graph, num_ev + 1);
	detail::Vec weights = detail::AssembleEdgeWeights(graph, solution);
	// create superpixel neighbourhood graph with edge strength
	SpectralGraph result(boost::num_vertices(graph));
	// FIXME proper edge indexing
	unsigned int eid_index = 0;
	for(auto eid : as_range(boost::edges(graph))) {
		auto edge = boost::add_edge(boost::source(eid, graph), boost::target(eid, graph), result);
		boost::put(boost::edge_weight, result, edge.first, weights[eid_index]);
		eid_index++;
	}
	return result;
}

}
