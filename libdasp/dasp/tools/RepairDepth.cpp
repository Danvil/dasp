/*
 * ReapairDepth.cpp
 *
 *  Created on: Mar 27, 2012
 *      Author: david
 */

#include "RepairDepth.hpp"
#include <Slimage/Parallel.h>
#include <iostream>

namespace dasp {

void RepairDepthFill(const slimage::Image1ui16& depth, const slimage::Image3ub& color)
{
	std::cout << "Repairing depth..." << std::endl;
	for(unsigned int i=0; i<depth.size(); i++) {
		if(depth[i] == 0) {
			depth[i] = 1000;
		}
	}

}

void RepairDepthGrow(const slimage::Image1ui16& depth, const slimage::Image3ub& color)
{
	slimage::Image1ui16 depth_raw = depth.clone();
	int w = depth.width();
	const unsigned int cSamples = 8;
	int neighbor_offsets[cSamples] = {
			-1-w, -w, +1-w,
			-1  ,     +1  ,
			-1+w, +w, +1+w
	};
	int weights[cSamples] = {
			1, 2, 1,
			2,    2,
			1, 2, 1
	};
	for(unsigned int y=1; y<depth.height()-1; y++) {
		for(unsigned int x=1; x<depth.width()-1; x++) {
			unsigned int center_index = x + y * w;
			if(depth[center_index] != 0) {
				continue;
			}
			int sum_w = 0;
			int sum_v = 0;
			for(unsigned int i=0; i<cSamples; i++) {
				uint16_t v = depth_raw[center_index + neighbor_offsets[i]];
				if(v == 0) continue;
				int w = weights[i];
				sum_w += w;
				sum_v += w * static_cast<int>(v);
			}
			depth(x,y) = (sum_w > 0) ? static_cast<uint16_t>(sum_v / sum_w) : 0;
		}
	}
}

void RepairDepth(const slimage::Image1ui16& depth, const slimage::Image3ub& color)
{
	RepairDepthGrow(depth, color);
	RepairDepthGrow(depth, color);
	RepairDepthGrow(depth, color);
}

void SmoothDepth(const slimage::Image1ui16& depth, const slimage::Image3ub& color)
{
	slimage::Image1ui16 depth_raw = depth.clone();
	int w = depth.width();
	const unsigned int cSamples = 9;
	int neighbor_offsets[cSamples] = {
			-1-w, -w, +1-w,
			-1  ,  0, +1  ,
			-1+w, +w, +1+w
	};
	int weights[cSamples] = {
			1, 2, 1,
			2, 4, 2,
			1, 2, 1
	};
	for(unsigned int y=1; y<depth.height()-1; y++) {
		for(unsigned int x=1; x<depth.width()-1; x++) {
			unsigned int center_index = x + y * w;
			int sum_w = 0;
			int sum_v = 0;
			for(unsigned int i=0; i<cSamples; i++) {
				uint16_t v = depth_raw[center_index + neighbor_offsets[i]];
				if(v == 0) continue;
				int w = weights[i];
				sum_w += w;
				sum_v += w * static_cast<int>(v);
			}
			depth(x,y) = (sum_w > 0) ? static_cast<uint16_t>(sum_v / sum_w) : 0;
		}
	}
}

}
