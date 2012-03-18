/*
 * AutoDepth.hpp
 *
 *  Created on: Feb 14, 2012
 *      Author: david
 */

#ifndef AUTODEPTH_HPP_
#define AUTODEPTH_HPP_
//----------------------------------------------------------------------------//
#include <ctype.h>
#include <vector>
//----------------------------------------------------------------------------//
namespace dasp {
//----------------------------------------------------------------------------//

namespace AutoFindDepthRange
{
	struct DepthHistogram
	{
		DepthHistogram(uint16_t depth_min, uint16_t depth_max, unsigned int bin_count)
		: depth_min_(depth_min), depth_max_(depth_max), bins_(bin_count, 0) {
			scl_ = float(bins_.size()) / float(depth_max_ - depth_min_);
		}

		void add(uint16_t x) {
			unsigned int bin;
			if(x < depth_min_) {
				bin = 0;
			}
			else if(x >= depth_max_) {
				bin = bins_.size() - 1;
			}
			else {
				float p = float(x - depth_min_) * scl_;
				bin = (unsigned int)p;
			}
			bins_[bin] ++;
		}

		unsigned int countBins() const {
			return bins_.size();
		}

		std::vector<unsigned int> getBins() const {
			return bins_;
		}

	private:
		uint16_t depth_min_, depth_max_;
		float scl_;
		std::vector<unsigned int> bins_;

	public:
		friend void Find(const DepthHistogram& h, uint16_t& min, uint16_t& max);
	};

	void Find(const DepthHistogram& h, uint16_t& min, uint16_t& max) {
//		std::vector<float> roi_depth_bin_cut(cRoiDepthBins);
//		unsigned int Ix = roi_depth_bins[0];
//		for(unsigned int i=1; i<cRoiDepthBins; i++) {
//			float t = float(i) / float(cRoiDepthBins - 1);
//			t *= t;
//			roi_depth_bin_cut[i] = float(Ix) / t;
//			Ix += roi_depth_bins[i];
//		}
//		auto cutit = std::max_element(roi_depth_bin_cut.begin(), roi_depth_bin_cut.end());
//		depth_max = ((cutit - roi_depth_bin_cut.begin()) * (cRoiDepthMax - cRoiDepthMin)) / cRoiDepthBins + cRoiDepthMin;

		unsigned int i = 1;
		while(h.bins_[i] < 10) i++;
		min = ((i-1) * (h.depth_max_ - h.depth_min_)) / h.countBins() + h.depth_min_;
		while(i+1 < h.countBins() && h.bins_[i] >= 10) i++;
		max = ((i+1) * (h.depth_max_ - h.depth_min_)) / h.countBins() + h.depth_min_;
	}
}

//----------------------------------------------------------------------------//
}
//----------------------------------------------------------------------------//
#endif
