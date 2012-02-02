/*
 * TreeReduction.cpp
 *
 *  Created on: Jan 29, 2012
 *      Author: david
 */


#include "TreeReduction.hpp"
#include <algorithm>
#include <map>
#include <cmath>

namespace Romeo {

namespace TreeReduction
{

	Graph MinimalCostEdges(const Graph& input)
	{
		unsigned int n = 3 * std::sqrt(input.edges.size());
		Graph result(input.nodes_);
		result.edges.resize(n);
		std::partial_sort_copy(input.edges.begin(), input.edges.end(), result.edges.begin(), result.edges.end(),
				[](const Edge& x, const Edge& y) { return x.cost < y.cost; });
		return result;
	}

	Graph MinimalSpanningCutting(const Graph& input, const float cut_param, std::vector<unsigned int>* labels)
	{
		struct Segment {
			float max_w;
			float tau;
			std::vector<unsigned int> vertices;
		};

		unsigned int nodes = input.nodes_;

		// initial segmentation each vertex on its own
		std::map<unsigned int,Segment> S;
		std::vector<unsigned int> seg_ids(nodes);
		for(unsigned int i=0; i<nodes; i++) {
			Segment& s = S[i];
			s.max_w = 0.0f;
			s.tau = cut_param;
			s.vertices = {i};
			seg_ids[i] = i;
		}

		// sort edges
		std::vector<Edge> edges = input.edges;
		std::sort(edges.begin(), edges.end(),
				[](const Edge& x, const Edge& y) { return x.cost < y.cost; });

		Graph graph(nodes);

		// iterate over all edges
		for(unsigned int k=0; k<edges.size(); k++) {
			const Edge& e = edges[k];
			unsigned int sid1 = seg_ids[e.a];
			unsigned int sid2 = seg_ids[e.b];
			// if vertices of edge in same segment -> use edge and continue
			if(sid1 == sid2) {
				graph.edges.push_back(e);
				continue;
			}
			// assert that sid1 < sid2 for minimal segment ids ... probably useless
			if(sid2 < sid1) {
				std::swap(sid1, sid2);
			}
			Segment& s1 = S[sid1];
			Segment& s2 = S[sid2];
			// test if we want to joint segments
			if(e.cost <= std::min(s1.max_w + s1.tau, s2.max_w + s2.tau)) {
				// join s2 into s1
				s1.max_w = std::max(e.cost, std::max(s1.max_w, s2.max_w));
				for(unsigned int vid : s2.vertices) {
					seg_ids[vid] = sid1;
				}
				s1.vertices.insert(s1.vertices.begin(), s2.vertices.begin(), s2.vertices.end());
				s1.tau = cut_param / s1.vertices.size();
				S.erase(sid2);
				graph.edges.push_back(e);
			}
		}

		if(labels) {
			*labels = seg_ids;
		}

		return graph;
	}

};

}
