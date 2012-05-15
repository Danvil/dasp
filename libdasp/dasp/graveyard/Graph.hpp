/*
 * TreeReduction.hpp
 *
 *  Created on: Jan 29, 2012
 *      Author: david
 */

#ifndef DASP_GRAPH_HPP_
#define DASP_GRAPH_HPP_

#include <boost/graph/edge_list.hpp>
#include <vector>
#include <functional>

namespace dasp {
namespace graph {

	struct Graph
	{
		struct Edge {
			float c_px;
			float c_world;
			float c_color;
			float c_normal;
			float cost;
		};

		typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, boost::no_property, Edge> graph_t;
		typedef graph_t::edge_descriptor edge_id_t;
		typedef graph_t::edge_iterator edge_it_t;

		Graph() { }

		Graph(unsigned int n) {
			for(unsigned int i=0; i<n; i++) {
				boost::add_vertex(graph_);
			}
		}

		unsigned int numNodes() const {
			return boost::num_vertices(graph_);
		}

		unsigned int numEdges() const {
			return boost::num_edges(graph_);
		}

		void add(unsigned int sid_a, unsigned int sid_b, const Edge& edge) {
			edge_id_t eid;
			bool ok;
			boost::tie(eid,ok) = boost::add_edge(i_keep, *it, graph_);
			graph_[eid] = edge;
		}

		void foreachEdge(const std::function<void(const Edge& edge)>& f) const {
			edge_it_t it, it_end;
			for(boost::tie(it, it_end)=boost::edges(graph_); it!=it_end; ++it) {
				return f(graph_[*it]);
			}
		}

		void mergeInto(unsigned int i_keep, unsigned int i_del) {
			graph_t::adjacency_iterator it, end;
			for(boost::tie(it,end)=boost::adjacent_vertices(i_del,graph_); it!=end; ++it) {
				boost::add_edge(i_keep, *it, graph_);
			}
			boost::clear_vertex(i_del, graph_);
			boost::remove_vertex(i_del, graph_);
		}

		friend Graph MinimalCostEdges(const Graph& input);

		friend Graph MinimalSpanningCutting(const Graph& input, const float cut_param, std::vector<unsigned int>* labels);

	private:
		void addImpl(const Edge& e) {
			edges_.push_back(e);
			edges_.back().id = edges_.size() - 1;
		}

	private:
		graph_t graph_;

	};

	Graph MinimalCostEdges(const Graph& input);

	Graph MinimalSpanningCutting(const Graph& input, const float cut_param, std::vector<unsigned int>* labels=0);

}}

#endif
