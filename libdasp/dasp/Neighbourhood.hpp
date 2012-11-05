/*
 * Neighbourhood.hpp
 *
 *  Created on: May 18, 2012
 *      Author: david
 */

#ifndef DASP_NEIGHBOURHOOD_HPP_
#define DASP_NEIGHBOURHOOD_HPP_

#include "Superpixels.hpp"
#include "Graph.hpp"
#include <Slimage/Slimage.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/copy.hpp>

namespace dasp
{

	/** Computes points which lie on the border between segments
	 * @param list of point indices
	 */
	std::vector<unsigned int> ComputeAllBorderPixels(const Superpixels& superpixels);

	struct NeighborGraphSettings
	{
		bool cut_by_spatial;
		float spatial_distance_mult_threshold;
		float pixel_distance_mult_threshold;
		float min_border_overlap;
		unsigned min_abs_border_overlap;

		static NeighborGraphSettings SpatialCut() {
			return NeighborGraphSettings{true, 5.0f, 4.0f, 0.0f, 1};
		}
		static NeighborGraphSettings NoCut() {
			return NeighborGraphSettings{false, 5.0f, 4.0f, 0.0f, 1};
		}
	};

	/** Computes the superpixel neighborhood graph
	 * Superpixels are neighbors if they share border pixels.
	 */
	BorderPixelGraph CreateNeighborhoodGraph(const Superpixels& superpixels, NeighborGraphSettings settings=NeighborGraphSettings::SpatialCut());

	/** Computes edge weights for a superpixel graph using the given metric
	 * Metric : Point x Point -> float (should be lightweight)
	 */
	template<typename Graph, typename Metric>
	EdgeWeightGraph ComputeEdgeWeights(const Superpixels& superpixels, const Graph& graph, const Metric& metric)
	{
		EdgeWeightGraph result;
		boost::copy_graph(graph, result, boost::edge_copy([&superpixels, &graph, &result, &metric](typename Graph::edge_descriptor src, EdgeWeightGraph::edge_descriptor dst) {
			unsigned int ea = source_superpixel_id(src, graph);
			unsigned int eb = target_superpixel_id(src, graph);
			float w = metric(superpixels.cluster[ea].center, superpixels.cluster[eb].center);
			boost::put(boost::edge_weight, result, dst, w);
		}));
		return result;
	}

	/** Computes a border image
	 * All border pixels for each edge are set to the edge weight in the image.
	 * Default value for non border pixels is 0.
	 * Border pixels are given in index form x + y*w where w must be the same given as parameter.
	 * @param w width of result image
	 * @param h width of result image
	 * @param graph a boost graph
	 * @param weights an weight property map of type float
	 * @param border_pixels an edge property map of type std::vector<unsigned int>
	 * @param return image with painted border pixels
	 */
	template<typename Graph, typename WeightMap, typename BorderPixelMap>
	slimage::Image1f CreateBorderPixelImage(unsigned int w, unsigned int h, const Graph& graph, WeightMap weights, BorderPixelMap border_pixels)
	{
		slimage::Image1f result(w, h, slimage::Pixel1f{0.0f});
		for(auto eid : as_range(boost::edges(graph))) {
			float v = boost::get(weights, eid);
			for(unsigned int pid : boost::get(border_pixels, eid)) {
				result[pid] = v;
			}
		}
		return result;
	}

	/** Smoothes a contour image and converts to unsigned char */
	slimage::Image1ub CreateSmoothedContourImage(const slimage::Image1f& src, float scl);

}

#endif
