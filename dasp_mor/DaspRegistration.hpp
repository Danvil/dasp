/*
 * DaspRegistration.hpp
 *
 *  Created on: Apr 17, 2012
 *      Author: david
 */

#ifndef DASPREGISTRATION_HPP_
#define DASPREGISTRATION_HPP_

#include "PointSet.hpp"
#include "ObjectModel.hpp"
#include <dasp/Superpixels.hpp>
#include <dasp/Plots.hpp>
#include <dasp/Segmentation.hpp>
#include <Slimage/Slimage.hpp>
#include <Slimage/IO.hpp>
#include <Slimage/Gui.hpp>
#include <Eigen/Dense>
#include <boost/bind.hpp>
#include <boost/format.hpp>
#include <boost/math/constants/constants.hpp>
#include <iostream>
#include <vector>

static constexpr bool cGui = true;


/** An Rgbd frame with color and depth image */
struct Rgbd
{
	slimage::Image3ub color;
	slimage::Image1ui16 depth;

	void load(const std::string& fn) {
		color = slimage::Load3ub(fn + "_color.png");
		depth = slimage::Load1ui16(fn + "_depth.pgm");
	}

	static Rgbd Load(const std::string& fn) {
		Rgbd x;
		x.load(fn);
		return x;
	}
};

/** A set of depth-adaptive superpixels */
struct DaspPointSet
{
	dasp::Superpixels superpixel;

	DaspPointSet() {
	}

	DaspPointSet(const Rgbd& u) {
		createWithDasp(u);
	}

	void createWithDasp(const Rgbd& u) {
		dasp::Parameters opt;
		opt.base_radius = 0.035f;
		opt.camera = dasp::Camera{320.0f, 240.0f, 540.0f, 0.001f};
		opt.is_repair_depth = true;
		superpixel = dasp::ComputeSuperpixels(u.color, u.depth, opt);
		// extract points
	}

	PointSet createPointSet() const {
		PointSet points(superpixel.clusterCount());
		for(unsigned int i=0; i<points.size(); i++) {
			points[i].position = superpixel.cluster[i].center.world;
			points[i].normal = superpixel.cluster[i].center.computeNormal();
			points[i].color = superpixel.cluster[i].center.color;
		}
		return points;
	}

};

constexpr int NoId = -1;

inline bool IsValidId(int id) {
	return id != NoId;
}

struct Partition
{
	std::vector<int> partition_uid_;

	unsigned int partition_count;

	int operator()(unsigned int i) const {
		return partition_uid_[i];
	}
};

struct Pairings
{
	struct Pair {
		unsigned int source_id;
		unsigned int target_id;
		float weight;
	};
	std::vector<std::vector<Pair>> pairings_;
};

struct Transformation
{
	std::vector<Eigen::Affine3f> T;

	const Eigen::Affine3f& operator[](std::size_t i) const {
		return T[i];
	}

	Eigen::Affine3f& operator[](std::size_t i) {
		return T[i];
	}

	PointSet transform(const PointSet& points, const Partition& p) const {
		PointSet q = points;
		for(unsigned int i=0; i<q.size(); i++) {
			int id = p(i);
			if(IsValidId(id)) {
				q[i] = T[id] * q[i];
			}
		}
		return q;
	}

	static Transformation Identity(unsigned int n=1) {
		return Transformation{std::vector<Eigen::Affine3f>(n, Eigen::Affine3f::Identity())};
	}

};

namespace impl {

	struct NearestPairing
	{
		Pairings operator()(const PointSet& pnts_source, const PointSet& pnts_target) {
			unsigned int n = pnts_source.size();
			std::vector<Pairings::Pair> u(n);
			for(unsigned int i=0; i<n; i++) {
				Pairings::Pair& p = u[i];
				p.source_id = i;
				p.target_id = FindClosestPoint(pnts_target, pnts_source[i]);
				p.weight = 1.0f;
			}
			return Pairings{{u}};
		}
	};

	struct PartitionNearestPairing
	{
		Pairings operator()(const PointSet& pnts_source, const PointSet& pnts_target, const Partition& partition) {
			unsigned int n = pnts_source.size();
			Pairings pairings;
			pairings.pairings_.resize(partition.partition_count);
			for(unsigned int i=0; i<n; i++) {
				Pairings::Pair p;
				p.source_id = i;
				p.target_id = FindClosestPoint(pnts_target, pnts_source[i]);
				p.weight = 1.0f;
				pairings.pairings_[partition(i)].push_back(p);
			}
			return pairings;
		}
	};
}

/** Performs iterative closest points for depth-adative superpixels */
struct IterativeClosestPoints
{
public:
	typedef boost::function<Pairings(const PointSet& pnts)> PairingFunctionType;

	IterativeClosestPoints() {
		is_ready = false;
	}

	const Transformation& getTransformation() const {
		return T;
	}

	const Pairings& getPairings() const {
		return pairings_;
	}

	const PointSet& getSourcePoints() const {
		return pnts_source_0;
	}

	const PointSet& getTargetPoints() const {
		return pnts_target;
	}

	void start(const PointSet& pnts_src, const PointSet& pnts_dst, const Transformation& T_initial, const Partition& partition, PairingFunctionType f) {
		pairing_function_ = f;
		pnts_source_0 = pnts_src;
		pnts_target = pnts_dst;
		T = T_initial;
		partition_ = partition;

		PointSet pnts_source = T.transform(pnts_source_0, partition_);
		pairings_ = pairing_function_(pnts_source);
	}

	bool isReady() const {
		return is_ready;
	}

	void step() {
		if(is_ready) {
			return;
		}
		PointSet pnts_source = T.transform(pnts_source_0, partition_);

		pairings_ = pairing_function_(pnts_source);

		unsigned int partition_count = pairings_.pairings_.size();
		std::vector<bool> ready(partition_count, false);
		for(unsigned int i=0; i<partition_count; i++) {
			if(ready[i]) {
				continue;
			}
			Eigen::Affine3f delta = IcpStepImpl(pnts_source, pnts_target, pairings_.pairings_[i]);
			T[i] = delta * T[i];
			if(IsNearIdentity(delta)) {
				std::cout << "Segment " << i << " converged" << std::endl;
				ready[i] = true;
			}
		}
		if(std::count_if(ready.begin(), ready.end(), [](bool q) { return q; }) == partition_count) {
			is_ready = true;
		}
	}

	void run() {
		while(!is_ready) {
			step();
		}
	}

private:
	PointSet pnts_source_0;
	PointSet pnts_target;
	PairingFunctionType pairing_function_;
	Partition partition_;
	Pairings pairings_;
	Transformation T;
	bool is_ready;

public:
	static bool IsNearIdentity(const Eigen::Affine3f& T) {
		const float cMaxDt = 0.0005f; // 0.5 mm
		const float cMaxDr = 0.5f / 180.0f * boost::math::constants::pi<float>(); // 0.5 deg
		float dt = T.translation().norm();
		float dr = Eigen::Quaternionf(T.rotation()).angularDistance(Eigen::Quaternionf::Identity());
		return dt < cMaxDt && dr < cMaxDr;
	}

	static std::vector<int> ComputePointCorrespondence(const PointSet& pnts_source, const PointSet& pnts_target)
	{
		unsigned int n = pnts_source.size();
		std::vector<int> cp_map(n);
		for(unsigned int i=0; i<n; i++) {
			cp_map[i] = FindClosestPoint(pnts_target, pnts_source[i]);
		}
		return cp_map;
	}

	/** Performs a step in iterative closest points
	 * @param T initial transformation
	 * @param pnts_source point set to move around
	 * @param pnts_target point set to converge to
	 * @param correspondence for each point in pnts_source the corresponding point in pnts_target
	 */
	static Eigen::Affine3f IcpStepImpl(const PointSet& pnts_source, const PointSet& pnts_target, const std::vector<Pairings::Pair>& pairings);

};

struct Frame
{
	Rgbd rgdb;

	DaspPointSet dasp;

	dasp::Segmentation segments;

	std::vector<slimage::Pixel3ub> segment_colors;

	Partition partition;

	Pairings pairings;

	Transformation T;

	unsigned int countSegments() const {
		return partition.partition_count;
	}

	unsigned int countSuperpixels() const {
		return dasp.superpixel.clusterCount();
	}

	PointSet transformPoints() const {
		return T.transform(dasp.createPointSet(), partition);
	}

	unsigned int countSegmentSuperpixels(unsigned int seg) const {
		return pairings.pairings_[seg].size();
	}

};

class IcpBatchObject;

struct DaspMultiIcp
{
	friend class IcpBatchObject;

public:
	DaspMultiIcp(const std::vector<Rgbd>& frames) {
		frames_.reserve(frames.size());
		current_frame_ = 0;
		for(Rgbd x : frames) {
			addFrame(x);
		}
	}

	void addFrame(const Rgbd& rgbd) {
		std::cout << "IcpBatchObject: Adding frame ..." << std::endl;
		Frame frame;
		frame.rgdb = rgbd;
		frame.dasp = DaspPointSet(rgbd);

		dasp::SpectralSettings segs_settings;
		segs_settings.num_eigenvectors = 24;
		segs_settings.w_spatial = 1.0f;
		segs_settings.w_color = 2.0f;
		segs_settings.w_normal = 3.0f;

		frame.segments = dasp::SpectralSegmentation(frame.dasp.superpixel, segs_settings);
		frame.segments.ucm(frame.dasp.superpixel, 5.0f); // FIXME what is the threshold?

		frame.partition.partition_count = frame.segments.segment_count;
		frame.partition.partition_uid_.resize(frame.segments.cluster_labels.size());
		std::copy(frame.segments.cluster_labels.begin(), frame.segments.cluster_labels.end(), frame.partition.partition_uid_.begin());

		frame.pairings.pairings_.resize(frame.partition.partition_count);
		for(unsigned int i=0; i<frame.countSuperpixels(); i++) {
			frame.pairings.pairings_[frame.partition(i)].push_back(Pairings::Pair{i,0,0.0f});
		}

//		// FIXME HACK HACK HACK for the case with only one partition
//		frame.partition = Partition{ std::vector<int>(dasp_points_.back().superpixel.clusterCount(), 0), 1};

		std::cout << "\tNumber of superpixels: " << frame.dasp.superpixel.clusterCount() << std::endl;
		std::cout << "\tNumber of segments: " << frame.partition.partition_count << std::endl;

		if(cGui) {
			slimage::gui::Show("color", frame.rgdb.color);
			slimage::gui::Show("dasp", dasp::plots::PlotClusters(frame.dasp.superpixel, dasp::plots::ClusterPoints, dasp::plots::Color));
			slimage::gui::Show("segments", frame.segments.computeLabelImage(frame.dasp.superpixel));
			slimage::gui::Show("contours", dasp::Segmentation::CreateSmoothedContourImage(frame.segments.boundaries, 0.013f));
			std::cout << "Press key to continue" << std::endl;
			slimage::gui::WaitForKeypress();
		}

		frame.T = Transformation::Identity(frame.partition.partition_count);

		//frame.segment_colors = micp_->dasp_segments_[micp_->current_frame_].computeSegmentColors(micp_->dasp_points_[micp_->current_frame_].superpixel);
		frame.segment_colors = dasp::plots::CreateRandomColors(frame.segments.segment_count);

		assert(frame.segment_colors.size() == frame.pairings.pairings_.size());

		frames_.push_back(frame);

		if(current_frame_ == 0 && world_.size() == 0) {
			updateWorld();
		}
	}

	void step() {
		if(icp_) {
			icp_->step();
			getFrame(current_frame_).T = icp_->getTransformation();
			getFrame(current_frame_).pairings = icp_->getPairings();
			if(icp_->isReady()) {
				updateWorld();
				writeGraphML("register.graphml");
				icp_.reset();
			}
		}
	}

	bool isFinished() const {
		return !icp_;
	}

	const Frame& getFrame(unsigned int i) const {
		return frames_[i];
	}

	Frame& getFrame(unsigned int i) {
		return frames_[i];
	}

	void writeGraphML(const std::string& fn);

private:
	void next() {
		if(frames_.size() < 2 || current_frame_ + 1 >= frames_.size()) {
			icp_.reset();
			return;
		}

		current_frame_ ++;

		icp_.reset(new IterativeClosestPoints());

		PointSet pnts_source = getFrame(current_frame_).dasp.createPointSet();
		PointSet pnts_target = getFrame(current_frame_ - 1).dasp.createPointSet();

		Transformation Tinital = Transformation::Identity(getFrame(current_frame_).countSegments());

//		impl::NearestPairing c_algo_;
//		IterativeClosestPoints::PairingFunctionType pairing_fnc = boost::bind(&impl::NearestPairing::operator(), c_algo_, _1, pnts_target);
		impl::PartitionNearestPairing c_algo_;
		IterativeClosestPoints::PairingFunctionType pairing_fnc = boost::bind(&impl::PartitionNearestPairing::operator(), c_algo_, _1, pnts_target, getFrame(current_frame_).partition);

		icp_->start(pnts_source, pnts_target, Tinital, getFrame(current_frame_).partition, pairing_fnc);
		std::cout << "IcpBatchObject: current frame=" << current_frame_ << ", number of points=" << pnts_source.size() << std::endl;
	}

	void updateWorld();

private:
	std::vector<Frame> frames_;

	unsigned int current_frame_;

	boost::shared_ptr<IterativeClosestPoints> icp_;

	ObjectModelGroup world_;

};

#endif
