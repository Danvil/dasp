/*
 * DaspRegistration.cpp
 *
 *  Created on: Apr 20, 2012
 *      Author: david
 */

#include "DaspRegistration.hpp"
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphml.hpp>
#include <boost/format.hpp>
#include <iostream>
#include <map>
#include <set>

Eigen::Affine3f IterativeClosestPoints::IcpStepImpl(const PointSet& pnts_source, const PointSet& pnts_target, const std::vector<Pairings::Pair>& pairings)
{
//	std::cout << "Pairings in segment: " << pairings.size() << std::endl;
	if(pairings.size() < 5) {
		return Eigen::Affine3f::Identity();
	}
	// assemble problem accordingly to
	// http://www.cs.princeton.edu/~smr/papers/icpstability.pdf
	typedef Eigen::Matrix<float,6,6> Mat6;
	typedef Eigen::Matrix<float,6,1> Vec6;
	Mat6 A = Mat6::Zero();
	Vec6 b = Vec6::Zero();
	for(const Pairings::Pair& p : pairings) {
		const Point& pa = pnts_source[p.source_id];
		const Point& pb = pnts_target[p.target_id];
		Eigen::Vector3f p = pa.position;
		Eigen::Vector3f q = pb.position;
		Eigen::Vector3f n = pa.normal;
		Eigen::Vector3f c = p.cross(n);
		Eigen::Vector3f d = p - q;
		float dn = d.dot(n);
		b[0] += c[0] * dn;
		b[1] += c[1] * dn;
		b[2] += c[2] * dn;
		b[3] += n[0] * dn;
		b[4] += n[1] * dn;
		b[5] += n[2] * dn;
		for(int u=0; u<3; u++) {
			for(int v=u; v<3; v++) {
				float x = c[u] * c[v];
				A(u,v) += x;
				A(v,u) += x;
			}
		}
		for(int u=0; u<3; u++) {
			for(int v=0; v<3; v++) {
				float x = c[u] * n[v];
				A(u,3+v) += x;
				A(3+v,u) += x;
			}
		}
		for(int u=0; u<3; u++) {
			for(int v=u; v<3; v++) {
				float x = n[u] * n[v];
				A(u+3,v+3) += x;
				A(v+3,u+3) += x;
			}
		}
	}
	// solve equations
//	std::cout << A << std::endl;
//	std::cout << b.transpose() << std::endl;
	Vec6 x = A.ldlt().solve(-b);
	// construct transformation matrix
	Eigen::Matrix3f Trot;
	Trot <<
			1, -x[2], x[1],
			x[2], 1, -x[0],
			-x[1], x[0], 1;
	Eigen::Vector3f Tpos;
	Tpos << x[3], x[4], x[5];
	Eigen::Affine3f Tinc;
	Tinc.linear() = Trot;
	Tinc.translation() = Tpos;
	return Tinc;
}

std::string ColorToHex(const slimage::Pixel3ub& color) {
	static boost::format hex("#%02x%02x%02x");
	return (hex % static_cast<int>(color[0]) % static_cast<int>(color[1]) % static_cast<int>(color[2])).str();
}

enum vertex_size_t { vertex_size };

namespace boost
{
    BOOST_INSTALL_PROPERTY(vertex, size);
}

void DaspMultiIcp::writeGraphML(const std::string& fn)
{
	using namespace boost;
	using namespace std;
	typedef adjacency_list< vecS, vecS, directedS,
			property<vertex_name_t,string, property<vertex_color_t,string, property<vertex_size_t,unsigned int> > >,
			property<edge_weight_t,unsigned int>
			> Graph;

	// connect all superpixels with via pairings
	Graph g;

	boost::format node_name("node_%02d_%04d");
	boost::format seg_name("seg_%02d_%03d");

	// add nodes
//	std::vector<std::vector<Graph::vertex_descriptor> > vids(frames_.size());
	std::vector<std::vector<Graph::vertex_descriptor> > seg_vids(frames_.size());
	for(unsigned int k=0; k<=current_frame_; k++) {
		// add segment nodes
		for(unsigned int i=0; i<frames_[k].countSegments(); i++) {
			Graph::vertex_descriptor vid = add_vertex(g);
			seg_vids[k].push_back(vid);
			// TODO add vertex properties
			put(vertex_name_t(), g, vid, (seg_name % k % i).str());
			put(vertex_color_t(), g, vid, ColorToHex(frames_[k].segment_colors[i]));
			put(vertex_size_t(), g, vid, frames_[k].countSegmentSuperpixels(i));
		}
//		// add superpixel nodes
//		vids[k].resize(frames_[k].countSuperpixels());
//		for(unsigned int i=0; i<vids[k].size(); i++) {
//			Graph::vertex_descriptor vid = add_vertex(g);
//			vids[k][i] = vid;
//			// TODO add vertex properties
//			put(vertex_name_t(), g, vid, (node_name % k % i).str());
//			put(vertex_color_t(), g, vid, ColorToHex(frames_[k].segment_colors[frames_[k].partition(i)]));
//		}
	}

//	// add segment to superpixel edges
//	for(unsigned int k=0; k<=current_frame_; k++) {
//		for(unsigned int i=0; i<vids[k].size(); i++) {
//			auto r = add_edge(vids[k][i], seg_vids[k][frames_[k].partition(i)], g);
//			assert(r.second);
//			Graph::edge_descriptor eid = r.first;
//			// TODO add edge properties
//			put(edge_weight_t(), g, eid, 0.0f);
//		}
//	}

//	// add superpixel frame to frame edges
//	for(unsigned int k=1; k<=current_frame_; k++) {
//		for(const std::vector<Pairings::Pair>& u : frames_[k].pairings.pairings_) {
//			for(const Pairings::Pair& p : u) {
//				auto r = add_edge(vids[k][p.source_id], vids[k-1][p.target_id], g);
//				assert(r.second);
//				Graph::edge_descriptor eid = r.first;
//				// TODO add edge properties
//				put(edge_weight_t(), g, eid, 0.0f);
//			}
//		}
//	}

	// add segment frame to frame edges
	for(unsigned int k=1; k<=current_frame_; k++) {
		std::vector<std::map<unsigned int,unsigned int>> refs(frames_[k].countSegments());
		for(const std::vector<Pairings::Pair>& u : frames_[k].pairings.pairings_) {
			for(const Pairings::Pair& p : u) {
				unsigned int s = frames_[k].partition(p.source_id);
				unsigned int t = frames_[k-1].partition(p.target_id);
				refs[s][t] ++;
			}
		}
		for(unsigned int i=0; i<refs.size(); i++) {
			for(std::pair<unsigned int, unsigned int> p : refs[i]) {
				auto r = add_edge(seg_vids[k][i], seg_vids[k-1][p.first], g);
				//double q = static_cast<double>(p.second) / static_cast<double>(frames_[k].countSegmentSuperpixels(i));
				//assert(q <= 1.0);
				put(edge_weight_t(), g, r.first, p.second);
			}
		}
	}

	std::cout << "Writing graph with " << num_vertices(g) << " vertices and " << num_edges(g) << " edges to file '" << fn << "'" << std::endl;

	// http://thirld.com/blog/2012/01/31/making-yed-import-labels-from-graphml-files/
	dynamic_properties dp;
	dp.property("text", get(vertex_name_t(), g));
	dp.property("color", get(vertex_color_t(), g));
	dp.property("size", get(vertex_size_t(), g));
	dp.property("weight", get(edge_weight_t(), g));

	std::ofstream ofs(fn);
	write_graphml(ofs, g, dp, true);
}

void DaspMultiIcp::updateWorld()
{
	std::cout << "updating world " << current_frame_ << std::endl;
	// points from the current frame where registered to match points from the previous frame

	const Frame& current_frame = frames_[current_frame_];
	// create new world based on current segments
	PointSet current_points = current_frame.dasp.createPointSet();
	ObjectModelGroup current_world(current_frame.countSegments());
	for(unsigned int i=0; i<current_points.size(); i++) {
		current_world[current_frame.partition(i)].add(current_points[i]);
	}

	// for the first frame we are ready
	if(current_frame_ == 0) {
		world_ = current_world;
		return;
	}

	ObjectModelGroup new_world(current_frame.countSegments());

	// split old segments
	// world_ contains the world model built from all previous frames
	// objects in the world represent segments from exactlye the previous step
	for(unsigned int prev_seg_id=0; prev_seg_id<world_.size(); prev_seg_id++) {
		const Frame& previous_frame = frames_[current_frame_-1];
		// find all new segments which want to get parts from me
		std::set<unsigned int> split_to_ids_set;
		for(const std::vector<Pairings::Pair>& u : current_frame.pairings.pairings_) {
			for(const Pairings::Pair& p : u) {
				// source = current frame
				// target = previous frame
				// if partition of target is myself, then partition of source wants a part of me
				if(previous_frame.partition(p.target_id) == prev_seg_id) {
					split_to_ids_set.insert(current_frame.partition(p.source_id));
				}
			}
		}
		// no segments ... -> skip and keep as is
		if(split_to_ids_set.size() == 0) {
			continue;
		}

		// prepare splitting
		// for each segment of the current frame which wants a part of this segment from the previous step
		// transform the points to from current to previous coordinate frame
		// so all computation is in wrt the previous step
		std::vector<unsigned int> split_to_ids(split_to_ids_set.begin(), split_to_ids_set.end());
		std::vector<PointSet> split_to_pnts(split_to_ids.size()); // in coordinates of previous frame
		for(unsigned int j=0; j<split_to_ids.size(); j++) {
			unsigned int id = split_to_ids[j];
			// transform segments to coordinate system of previous step
			split_to_pnts[j] = current_frame.T[id] * current_world[id].points_;
		}

		// split
		ObjectModelGroup parts = world_[prev_seg_id].split(split_to_pnts);
		assert(parts.size() == split_to_ids.size());
		assert(parts.size() == split_to_pnts.size());
		// assign to new world
		for(unsigned int j=0; j<split_to_ids.size(); j++) {
			unsigned int nid = split_to_ids[j];
//			std::cout << "\tsplitting to " << nid << std::endl;
			new_world[nid].add(parts[j].points_);
		}
	}
	// transform segments to current frame
	world_.clear();
	world_.resize(current_world.size());
	for(unsigned int i=0; i<world_.size(); i++) {
		world_[i].points_ = current_frame.T[i].inverse() * new_world[i].points_;
//		world_[i].add(current_world[i].points_);
	}
}

