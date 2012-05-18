/*
 * Graph.hpp
 *
 *  Created on: May 18, 2012
 *      Author: david
 */

#ifndef DASP_GRAPH_HPP_
#define DASP_GRAPH_HPP_

#include <boost/graph/adjacency_list.hpp>

namespace dasp
{
	namespace detail
	{
		template<class Iter>
		struct iter_pair_range
		: std::pair<Iter,Iter>
		{
			iter_pair_range(const std::pair<Iter,Iter>& x)
			: std::pair<Iter,Iter>(x)
			{}

			Iter begin() const {
				return this->first;
			}

			Iter end() const {
				return this->second;
			}
		};
	}

	template<class Iter>
	inline detail::iter_pair_range<Iter> as_range(const std::pair<Iter,Iter>& x) {
		return detail::iter_pair_range<Iter>(x);
	}

	template<class Graph>
	inline detail::iter_pair_range<typename Graph::edge_iterator> edges_range(const Graph& graph) {
		return as_range(boost::edges(graph));
	}

	template<class Graph>
	inline detail::iter_pair_range<typename Graph::vertex_iterator> vertices_range(const Graph& graph) {
		return as_range(boost::vertices(graph));
	}

	struct borderpixels_t {
		typedef boost::edge_property_tag kind;
	};

	struct superpixel_id_t {
		typedef boost::vertex_property_tag kind;
	};

	typedef unsigned int SuperpixelId;

	typedef boost::property<boost::edge_weight_t, float> EdgeWeightProperty;

	typedef boost::property<borderpixels_t, std::vector<unsigned int>> EdgeBorderPixelsProperty;

	typedef boost::property<superpixel_id_t, SuperpixelId> VertexSuperpixelIdProperty;

	/** An undirected graph with a list of image border pixels per edge */
	typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, VertexSuperpixelIdProperty, EdgeBorderPixelsProperty> BorderPixelGraph;

	/** An undirected weighted graph */
	typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, VertexSuperpixelIdProperty, EdgeWeightProperty> EdgeWeightGraph;

	/** Gets the superpixel id of the source vertex of a superpixel graph edge
	 * SuperpixelGraph must have a superpixel_id_t edge property.
	 */
	template<typename SuperpixelGraph>
	inline SuperpixelId source_superpixel_id(const typename SuperpixelGraph::edge_descriptor& eid, const SuperpixelGraph& graph) {
		return boost::get(superpixel_id_t(), graph, boost::source(eid, graph));
	}

	/** Gets the superpixel id of the target vertex of a superpixel graph edge
	 * SuperpixelGraph must have a superpixel_id_t edge property.
	 */
	template<typename SuperpixelGraph>
	inline SuperpixelId target_superpixel_id(const typename SuperpixelGraph::edge_descriptor& eid, const SuperpixelGraph& graph) {
		return boost::get(superpixel_id_t(), graph, boost::source(eid, graph));
	}

	template<typename SuperpixelGraph>
	inline void put_superpixel_id(SuperpixelGraph& graph, const typename SuperpixelGraph::vertex_descriptor& vid, SuperpixelId sid) {
		boost::put(superpixel_id_t(), graph, vid, sid);
	}

	template<typename SuperpixelGraph>
	inline SuperpixelId get_superpixel_id(const SuperpixelGraph& graph, const typename SuperpixelGraph::vertex_descriptor& vid) {
		return boost::get(superpixel_id_t(), graph, vid);
	}

	namespace detail
	{
		template<typename SuperpixelGraph>
		SuperpixelGraph CreateSuperpixelGraph(unsigned int num_vertices) {
			SuperpixelGraph graph(num_vertices);
			for(auto vid : as_range(boost::vertices(graph))) {
				put_superpixel_id(graph, vid, static_cast<SuperpixelId>(vid));
			}
//			// more correct but slower
//			Graph g;
//			for(Superpixels i=0; i<superpixels.clusterCount(); i++) {
//				auto vid = boost::add_vertex(g);
//				boost::put(superpixelid_t(), g, vid, i);
//			}
			return graph;
		}
	}
}

#endif
