/*
 * Neighbourhood.cpp
 *
 *  Created on: May 18, 2012
 *      Author: david
 */

#include "Neighbourhood.hpp"
#include <Slimage/Convert.hpp>

namespace dasp
{

//		/** Computes a list of all pixels which are have label cid and have a face neighbor which has label cjd */
//		std::vector<unsigned int> ComputeBorderPixelsImpl(unsigned int cid, unsigned int cjd, const Superpixels& spc, const slimage::Image1i& labels);

//		template<typename Graph>
//		std::vector<std::vector<unsigned int> > ComputeBorderPixels(const Superpixels& superpixels, const Graph& graph) const {
//			slimage::Image1i labels = superpixels.ComputeLabels();
//			std::vector<std::vector<unsigned int> > borders;
//			borders.reserve(boost::num_edges(graph));
//			for(auto it=boost::edges(graph); it.first!=it.second; ++it.first) {
//				typename Graph::edge_descriptor eid = *it.first;
//				// compute pixels which are at the border between superpixels e.a and e.b
//				unsigned int i = boost::get( boost::source(eid, graph);
//				unsigned int j = boost::target(eid, graph);
//				// superpixel i should have less points than superpixel j
//				if(superpixels.cluster[i].pixel_ids.size() > superpixels.cluster[j].pixel_ids.size()) {
//					std::swap(i,j);
//				}
//				// find border pixels
//				borders.push_back( detail::ComputeBorderPixelsImpl(i, j, *this, labels) );
//			}
//			return borders;
//		}

std::vector<unsigned int> ComputeBorderPixelsImpl(unsigned int cid, unsigned int cjd, const Superpixels& spc, const slimage::Image1i& labels)
{
	const int w = static_cast<int>(spc.width());
	const int h = static_cast<int>(spc.height());
	const int d[4] = { -1, +1, -w, +w };
	std::vector<unsigned int> border;
	for(unsigned int pid : spc.cluster[cid].pixel_ids) {
		int x = pid % w;
		int y = pid / w;
		if(1 <= x && x+1 <= w && 1 <= y && y+1 <= h) {
			for(int i=0; i<4; i++) {
				int label = labels[pid + d[i]];
				if(label == static_cast<int>(cjd)) {
					border.push_back(pid);
				}
			}
		}
	}
	return border;
}

std::vector<int> ComputeBorderLabels(unsigned int cid, const Superpixels& spc, const slimage::Image1i& labels)
{
	const int w = static_cast<int>(spc.width());
	const int h = static_cast<int>(spc.height());
	const int d[4] = { -1, +1, -w, +w };
	std::vector<int> border;
	for(unsigned int pid : spc.cluster[cid].pixel_ids) {
		int x = pid % w;
		int y = pid / w;
		if(1 <= x && x+1 <= w && 1 <= y && y+1 <= h) {
			for(int i=0; i<4; i++) {
				int label = labels[pid + d[i]];
				if(label != cid && label != -1) {
					border.push_back(label);
				}
			}
		}
	}
	return border;
}

std::vector<std::vector<int>> ComputeBorders(const Superpixels& spc)
{
	slimage::Image1i labels = spc.ComputeLabels();
	std::vector<std::vector<int>> border_pixels(spc.cluster.size());
	for(unsigned int cid=0; cid<spc.cluster.size(); cid++) {
		border_pixels[cid] = ComputeBorderLabels(cid, spc, labels);
	}
	return border_pixels;
}

std::vector<unsigned int> ComputeAllBorderPixels(const Superpixels& superpixels)
{
	slimage::Image1i labels = superpixels.ComputeLabels();
	std::set<unsigned int> u;
	int face_neighbours[] = { -1, +1, -labels.width(), +labels.width() };
	for(unsigned int y=1; y<labels.height()-1; y++) {
		for(unsigned int x=1; x<labels.width()-1; x++) {
			int q = labels.index(x,y);
			unsigned int lc = labels[q];
			for(unsigned int i=0; i<4; i++) {
				if(labels[q + face_neighbours[i]] != lc) {
					u.insert(q);
				}
			}
		}
	}
	return std::vector<unsigned int>(u.begin(), u.end());
}

BorderPixelGraph CreateNeighborhoodGraph(const Superpixels& superpixels, NeighborGraphSettings settings)
{
	// create one node for each superpixel
	BorderPixelGraph neighbourhood_graph = detail::CreateSuperpixelGraph<BorderPixelGraph>(superpixels.clusterCount());
	// compute superpixel borders
	std::vector<std::vector<int> > borders = ComputeBorders(superpixels);
	// connect superpixels
	const float node_distance_threshold = settings.max_spatial_distance_mult * superpixels.opt.base_radius;
	for(unsigned int i=0; i<superpixels.cluster.size(); i++) {
		for(unsigned int j=i+1; j<superpixels.cluster.size(); j++) {
			if(settings.cut_by_spatial) {
				float d = (superpixels.cluster[i].center.world - superpixels.cluster[j].center.world).norm();
				// only test if distance is smaller than threshold
				if(d > node_distance_threshold) {
					continue;
				}
			}
			// compute intersection of border pixels
			std::vector<unsigned int> common_border(borders[i].size() + borders[j].size());
			auto common_border_end = std::set_intersection(borders[i].begin(), borders[i].end(), borders[j].begin(), borders[j].end(), common_border.begin());
			common_border.resize(common_border_end - common_border.begin());
			// test if superpixels have a common border
			unsigned int common_border_size = common_border.size();
			if(common_border_size < settings.min_abs_border_overlap) {
				continue;
			}
			float p = static_cast<float>(common_border_size) / static_cast<float>(std::min(borders[i].size(),borders[j].size()));
			if(p < settings.min_border_overlap) {
				continue;
			}
			// add edge
			BorderPixelGraph::edge_descriptor eid;
			bool ok;
			boost::tie(eid,ok) = boost::add_edge(i, j, neighbourhood_graph); // FIXME correctly convert superpixel_id to vertex descriptor
			boost::put(borderpixels_t(), neighbourhood_graph, eid, common_border);
//			NeighbourhoodGraphEdgeData& edge = neighbourhood_graph[eid];
//			// compute cost using color and normal
//			edge.c_px = metric::ImageDistanceRaw(cluster[i].center, cluster[j].center);
//			edge.c_world = metric::SpatialDistanceRaw(cluster[i].center, cluster[j].center) / (opt.base_radius * opt.base_radius); // F IXME HAAAACK
//			edge.c_color = metric::ColorDistanceRaw(cluster[i].center, cluster[j].center);
//			edge.c_normal = metric::NormalDistanceRaw(cluster[i].center, cluster[j].center);
//			// F IXME metric needs central place!
//			if(opt.weight_image == 0.0f) {
//				// dasp
//				MetricDASP fnc(opt.weight_spatial, opt.weight_color, opt.weight_normal, opt.base_radius);
//				edge.weight = fnc(cluster[i].center, cluster[j].center);
//			}
//			else {
//				// slic
//				MetricSLIC fnc(opt.weight_image, opt.weight_color);
//				edge.weight = fnc(cluster[i].center, cluster[j].center);
//			}
		}
	}
	return neighbourhood_graph;
}

slimage::Image1ub CreateSmoothedContourImage(const slimage::Image1f& src, float scl)
{
	slimage::Image1f tmp(src.dimensions(), slimage::Pixel1f{1.0f});
	for(unsigned int x=1; x<tmp.width()-1; x++) {
		for(unsigned int y=1; y<tmp.height()-1; y++) {
			float nb = src(x-1,y) + src(x+1,y) + src(x,y-1) + src(x,y+1);
			float v = src(x,y) + nb / 4.0f;
			tmp(x,y) = 1.0f - scl * v;
		}
	}
	slimage::Image1ub vis(tmp.dimensions());
	slimage::conversion::Convert(tmp, vis);
	return vis;
}

}
