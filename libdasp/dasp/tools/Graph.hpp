/*
 * TreeReduction.hpp
 *
 *  Created on: Jan 29, 2012
 *      Author: david
 */

#ifndef DASP_GRAPH_HPP_
#define DASP_GRAPH_HPP_

#include <vector>

namespace dasp {
namespace graph {

	struct Edge
	{
		unsigned int a, b;
		float c_world;
		float c_color;
		float c_normal;
		float cost;
	};

	struct Graph {
		Graph() : nodes_(0) {}

		Graph(unsigned int n) : nodes_(n) {}

		unsigned int nodes_;

		std::vector<Edge> edges;

		void mergeInto(unsigned int i_keep, unsigned int i_del) {
			for(unsigned int i=0; i<edges.size(); i++) {
				if(edges[i].a == i_del) {
					edges[i].a = i_keep;
				}
				if(edges[i].b == i_del) {
					edges[i].b = i_keep;
				}
			}
			removeDuplicatedEdges();
		}

		void removeDuplicatedEdges() {
			std::vector<Edge> neu;
			for(unsigned int i=0; i<edges.size(); i++) {
				bool keep = true;
				for(unsigned int j=i+1; j<edges.size(); j++) {
					if((edges[i].a == edges[j].a && edges[i].b == edges[j].b) || (edges[i].a == edges[j].b && edges[i].b == edges[j].a)) {
						keep = false;
						break;
					}
				}
				if(keep) {
					neu.push_back(edges[i]);
				}
			}
			edges = neu;
		}
	};

	Graph MinimalCostEdges(const Graph& input);

	Graph MinimalSpanningCutting(const Graph& input, const float cut_param, std::vector<unsigned int>* labels=0);

}}

#endif
