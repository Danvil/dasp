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

	struct Edge
	{
		unsigned int id;
		unsigned int a, b;
		float c_px;
		float c_world;
		float c_color;
		float c_normal;
		float cost;
	};

	struct Graph
	{
		Graph()
		: num_nodes_(0) {
		}

		Graph(unsigned int n)
		: num_nodes_(n) {
		}

		unsigned int numNodes() const {
			return num_nodes_;
		}

		unsigned int numEdges() const {
			return edges_.size();
		}

		const std::vector<Edge>& getEdges() const {
			return edges_;
		}

		std::vector<Edge>& getEdges() {
			return edges_;
		}

		void add(const Edge& edge) {
			addImpl(edge);
		}

		void foreachEdge(const std::function<void(const Edge& edge)>& f) const {
			for(const Edge& e : edges_) {
				return f(e);
			}
		}

		void mergeInto(unsigned int i_keep, unsigned int i_del) {
			for(unsigned int i=0; i<edges_.size(); i++) {
				if(edges_[i].a == i_del) {
					edges_[i].a = i_keep;
				}
				if(edges_[i].b == i_del) {
					edges_[i].b = i_keep;
				}
			}
			removeDuplicatedEdges();
		}

		void removeDuplicatedEdges() {
			std::vector<Edge> old = edges_;
			edges_.clear();
			edges_.reserve(old.size());
			for(unsigned int i=0; i<old.size(); i++) {
				bool keep = true;
				for(unsigned int j=i+1; j<old.size(); j++) {
					if((old[i].a == old[j].a && old[i].b == old[j].b) || (old[i].a == old[j].b && old[i].b == old[j].a)) {
						keep = false;
						break;
					}
				}
				if(keep) {
					addImpl(old[i]);
				}
			}
		}

		friend Graph MinimalCostEdges(const Graph& input);

		friend Graph MinimalSpanningCutting(const Graph& input, const float cut_param, std::vector<unsigned int>* labels);

	private:
		void addImpl(const Edge& e) {
			edges_.push_back(e);
			edges_.back().id = edges_.size() - 1;
		}

	private:
		unsigned int num_nodes_;

		std::vector<Edge> edges_;

	};

	Graph MinimalCostEdges(const Graph& input);

	Graph MinimalSpanningCutting(const Graph& input, const float cut_param, std::vector<unsigned int>* labels=0);

}}

#endif
