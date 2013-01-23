/*
 * SolveLapackTemplate.hpp
 *
 *  Created on: Jan 19, 2012
 *      Author: david
 */

#ifndef DASP_SPECTRAL_SOLVELAPACKTEMPLATE_HPP_
#define DASP_SPECTRAL_SOLVELAPACKTEMPLATE_HPP_

#include "../Common.hpp"
#include "../as_range.hpp"
#include <boost/graph/adjacency_list.hpp>
#include <iostream>
#include <fstream>
#include <vector>

extern "C" int ssyevr_(
		char*,char*,char*,
		int*,float*,int*,
		float*,float*,int*,int*,
		float*,
		int*,float*,float*,int*,int*,
		float*,int*,int*,int*,
		int*);

namespace graphseg { namespace detail {

template<typename Graph>
std::vector<EigenComponent> SolveLapackTemplate(const Graph& graph, unsigned int num_ev)
{

	unsigned int dim = boost::num_vertices(graph);
#ifdef SPECTRAL_VERBOSE
	unsigned int num_edges = boost::num_edges(graph);
	std::cout << "SpectralSegmentation: dimension=" << dim
			<< ", num_edges=" << num_edges
			<< ", matrix non-zero elements = " << 100*static_cast<float>(2 * num_edges) / static_cast<float>(dim*dim) << "%" << std::endl;
	std::vector<int> nodes_with_no_connection;
#endif
	// creating matrices
	Mat W = Mat::Zero(dim,dim);
	std::vector<float> Di(dim, 0.0f);
	for(auto eid : as_range(boost::edges(graph))) {
		unsigned int ea = boost::source(eid, graph);
		unsigned int eb = boost::target(eid, graph);
		float ew = boost::get(boost::edge_weight, graph, eid);
		if(std::isnan(ew)) {
			std::cerr << "Weight for edge (" << ea << "," << eb << ") is nan!" << std::endl;
			continue;
		}
		if(ew < 0) {
			std::cerr << "Weight for edge (" << ea << "," << eb << ") is negative!" << std::endl;
			continue;
		}
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
#ifdef SPECTRAL_VERBOSE
			nodes_with_no_connection.push_back(i);
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
#ifdef SPECTRAL_VERBOSE
	if(!nodes_with_no_connection.empty()) {
		std::cout << "Nodes without connections (#=" << nodes_with_no_connection.size() << "): ";
		for(int i : nodes_with_no_connection) {
			std::cout << i << ", ";
		}
		std::cout << std::endl;
	}
#endif	
	// compute matrix A = D^{-1/2} * (D - W) * D^{-1/2}
	Eigen::VectorXf di_sqrt_inv(dim);
	for(unsigned int i=0; i<dim; i++) {
		di_sqrt_inv[i] = 1.0f / std::sqrt(Di[i]);
	}
	Eigen::MatrixXf A(dim, dim);
	for(unsigned int i=0; i<dim; i++) {
		for(unsigned int j=0; j<dim; j++) {
			A(j,i) = - W(j,i) * di_sqrt_inv[i] * di_sqrt_inv[j];
		}
		A(i,i) = 1.0f;
	}

#ifdef SPECTRAL_VERBOSE
	{	// DEBUG
		std::ofstream ofs_A("/tmp/spectral_A.csv");
		for(unsigned int i=0; i<dim; i++) {
			for(unsigned int j=0; j<dim; j++) {
				ofs_A << A(i,j);
				if(j+1 == dim) {
					ofs_A << std::endl;
				}
				else {
					ofs_A << ", ";
				}
			}
		}
	}	// DEBUG
#endif

	// calling lapack!
	std::cout << "LAPACK: Start!" << std::endl;

	int N = dim;
	std::cout << "N=" << N << std::endl;
	float* data = A.data();
	float vl, vu;
	int il = 1, iu = 25;
	float accuracy = 0.00001f;
	int result_num_ew_found;
	float* result_ew = new float[N];
	float* result_ev = new float[N*N];
	int* result_isuppz = new int[2*N];
	int work_dim = -1;
	float* work = new float[1];
	int iwork_dim = -1;
	int* iwork = new int[1];
	int info;

	// ssyevr_(
	// 	"N", // JOBZ eigenvalues + eigenvectors
	// 	"A", // RANGE only some eigenvalues
	// 	"U", // UPLO upper triangle is stored
	// 	&N, // N order of A
	// 	data, // A upper triangle of A
	// 	&N, // LDA leading dimension of A
	// 	&vl, &vu, // VL,VU not used
	// 	&il, &iu, // IL, IU  range of eigenvalues returned
	// 	&accuracy, // ABSTOL accuracy
	// 	&result_num_ew_found, // M number of eigenvalues found
	// 	result_ew, // W computed eigenvalues
	// 	result_ev, // Z computed eigenvectors
	// 	&N, // LDZ leading dimension of Z
	// 	result_isuppz, // ISUPPZ
	// 	work, // WORK
	// 	&work_dim, // LWORK
	// 	iwork, // IWORK
	// 	&iwork_dim, // LIWORK
	// 	&info // INFO
	// 	);
	// work_dim = (int)*work;
	// iwork_dim = *iwork;

	work_dim = N*N;
	iwork_dim = N*N;

	delete[] work;
	delete[] iwork;
	std::cout << "LAPACK: work_dim (opt) =" << work_dim << std::endl;
	std::cout << "LAPACK: iwork_dim (opt) =" << iwork_dim << std::endl;
	work = new float[work_dim];
	iwork = new int[iwork_dim];

	ssyevr_(
		"N", // JOBZ eigenvalues + eigenvectors
		"A", // RANGE only some eigenvalues
		"U", // UPLO upper triangle is stored
		&N, // N order of A
		data, // A upper triangle of A
		&N, // LDA leading dimension of A
		0, 0, // VL,VU not used
		&il, &iu, // IL, IU  range of eigenvalues returned
		&accuracy, // ABSTOL accuracy
		&result_num_ew_found, // M number of eigenvalues found
		result_ew, // W computed eigenvalues
		result_ev, // Z computed eigenvectors
		&N, // LDZ leading dimension of Z
		result_isuppz, // ISUPPZ
		work, // WORK
		&work_dim, // LWORK
		iwork, // IWORK
		&iwork_dim, // LIWORK
		&info // INFO
		);

	std::cout << "LAPACK: Finished!" << std::endl;
	std::cout << "Number of eigenvalues: " << result_num_ew_found << std::endl;
	std::cout << "Info: " << info << std::endl;

	// return eigenvectors and eigenvalues
	std::vector<EigenComponent> solution(std::min<unsigned int>(num_ev, result_num_ew_found));
	for(std::size_t i=0; i<solution.size(); i++) {
		solution[i].eigenvalue = result_ew[i];
#ifdef SPECTRAL_VERBOSE
		std::cout << "SpectralSegmentation: eigenvalue #" << i << "=" << solution[i].eigenvalue << std::endl;
#endif
		solution[i].eigenvector = Eigen::VectorXf(dim);
		for(unsigned int j=0; j<dim; j++) {
			solution[i].eigenvector[j] = result_ev[i*dim + j]; // FIXME correct order?
		}
	}

	delete[] result_ew;
	delete[] result_ev;
	delete[] result_isuppz;
	delete[] work;
	delete[] iwork;

	return solution;
}

}}

#endif
