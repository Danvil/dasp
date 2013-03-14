#include "eval.hpp"
#include <dasp/Superpixels.hpp>
#include <boost/math/constants/constants.hpp>

namespace dasp {
namespace eval {

slimage::Image1i MarkBoundary(const slimage::Image1i& labels)
{
	slimage::Image1i boundary(labels.dimensions(), slimage::Pixel1i{1});
	for(unsigned int y=1; y+1<boundary.height(); y++) {
		for(unsigned int x=1; x+1<boundary.width(); x++) {
			int q = labels(x, y);
			int qp0 = labels(x+1, y  );
			int qm0 = labels(x-1, y  );
			int q0p = labels(x  , y+1);
			int q0m = labels(x  , y-1);
			if(q == qp0 && q == qm0 && q == q0p && q == q0m) {
				boundary(x,y) = 0;
			}
		}
	}
	return boundary;
}

std::pair<float,std::vector<float>> IsoperimetricQuotient(const Superpixels& u)
{
	// compute labels
	slimage::Image1i labels = u.ComputeLabels();
	// compute boundary image
	slimage::Image1i boundary = MarkBoundary(labels);
	// compute total pixel count and number of boundary pixels for each cluster
	std::vector<std::pair<int,int>> ipq_els(u.clusterCount());
	for(std::size_t i=0; i<labels.size(); i++) {
		int label = labels[i];
		if(label == -1) {
			continue;
		}
		ipq_els[label].first ++;	
		if(boundary[i] == 1) {
			ipq_els[label].second ++;
		}
	}
	// finalize
	// per cluster: v_i = 4*pi*A_i/L_i^2
	// total: sum_i v_i*A_i
	std::vector<float> ipq(ipq_els.size());
	float ipq_total = 0.0f;
	unsigned int num_total = 0;
	for(std::size_t i=0; i<ipq.size(); i++) {
		std::pair<int,int> q = ipq_els[i];
		float x = 4.0f*boost::math::constants::pi<float>()
			*static_cast<float>(q.first)
			/ static_cast<float>(q.second*q.second);
		ipq[i] = x;
		ipq_total += x*static_cast<float>(q.second);
		num_total += q.second;
	}
	return {ipq_total / static_cast<float>(num_total), ipq};
}

}}
