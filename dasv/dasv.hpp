#ifndef DASV_DASV_HPP
#define DASV_DASV_HPP

#include <rgbd.hpp>
#include <dasp/Array.hpp>
#include <dasp/Point.hpp>
#include <Slimage/Slimage.hpp>
#include <Eigen/Dense>
#include <boost/graph/adjacency_list.hpp>
#include <boost/iterator/transform_iterator.hpp>
#include <boost/function.hpp>
#include <vector>
#include <ostream>
#include <string>
#include <memory>

namespace dasv
{
	/** Sets a callback functions for displaying images */
	void DebugSetDisplayImageCallback(boost::function<void(const std::string& tag, const slimage::Image3ub& img)> f);

	/** Displays an image */
	void DebugDisplayImage(const std::string& tag, const slimage::Image3ub& img);

	/** Displays a matrix */
	void DebugDisplayImage(const std::string& tag, const Eigen::MatrixXf& mat, float min, float max);

	// /** Voxel data */
	typedef dasp::Point Point;

	/** A frame of voxels at a given timestamp */
	typedef dasp::ImagePoints RgbdData;

	typedef int cluster_id_type;

	constexpr cluster_id_type INVALID_CLUSTER_ID = -1;

	/** A voxel cluster (supervoxel) */
	struct Cluster
	{
		Eigen::Vector2f pixel;
		Eigen::Vector3f color;
		Eigen::Vector3f position;
		Eigen::Vector3f normal;
		float cluster_radius_px;
		int time;
		cluster_id_type cluster_id;
		bool valid;
		int label;

		friend std::ostream& operator<<(std::ostream& os, const Cluster& c)
		{
			os << c.time
				<< "\t" << c.cluster_id
				<< "\t" << (c.valid ? 1 : 0)
				<< "\t" << c.label
				<< "\t" << c.cluster_radius_px
				<< "\t" << c.pixel.x() << "\t" << c.pixel.y()
				<< "\t" << c.color.x() << "\t" << c.color.y() << "\t" << c.color.z()
				<< "\t" << c.position.x() << "\t" << c.position.y() << "\t" << c.position.z()
				<< "\t" << c.normal.x() << "\t" << c.normal.y() << "\t" << c.normal.z();
			return os;
		}

	};

	struct ClusterContainer
	{
		const std::vector<Cluster>& data() const { return data_; }
		std::size_t size() const { return data_.size(); }
		const Cluster& at(cluster_id_type i) const { return data_[i]; }
		Cluster& at(cluster_id_type i) { return data_[i]; }
		cluster_id_type add(const Cluster& c) {
			data_.push_back(c);
			cluster_id_type& id = data_.back().cluster_id;
			id = data_.size() - 1;
			return id;
		}
		cluster_id_type add() { return add(Cluster()); }
	private:
		std::vector<Cluster> data_;
	};

	struct ClusterList
	{
	public:
		static std::shared_ptr<ClusterContainer> s_storage_;
		typedef std::vector<cluster_id_type> index_container_t;
		struct Access {
			typedef Cluster& result_type;
			Cluster& operator()(cluster_id_type i) const { return s_storage_->at(i); }
		};
		typedef boost::transform_iterator<Access,index_container_t::const_iterator> cit_t;
		typedef boost::transform_iterator<Access,index_container_t::iterator> it_t;
		const index_container_t&  indices() const { return indices_; }
		index_container_t&  indices() { return indices_; }
		cit_t begin() const { return cit_t(indices_.cbegin(), Access()); }
		it_t begin() { return it_t(indices_.begin(), Access()); }
		cit_t end() const { return cit_t(indices_.cend(), Access()); }
		it_t end() { return it_t(indices_.end(), Access()); }
		std::size_t size() const { return indices_.size(); }
		Cluster& addCluster(const Cluster& c) {
			auto cid = s_storage_->add(c);
			indices_.push_back(cid);
			return s_storage_->at(cid);
		}
		Cluster& addCluster() { return addCluster(Cluster()); }
		Cluster& operator[](int i) { return s_storage_->at(indices_[i]); }
		const Cluster& operator[](int i) const { return s_storage_->at(indices_[i]); }
	private:
		index_container_t indices_;
	};

	/** Pixel to cluster assignment */
	struct Assignment
	{
		cluster_id_type cluster_id;
		float distance;

		bool hasValidCluster() const {
			return cluster_id != INVALID_CLUSTER_ID;
		}

		const Cluster& getCluster() const {
			return ClusterList::s_storage_->at(cluster_id);
		}

		Cluster& getCluster() {
			return ClusterList::s_storage_->at(cluster_id);
		}

		static Assignment Empty() {
			return { INVALID_CLUSTER_ID, 10000000.0f };
		}
	};

	/** Frame pixel to cluster assignment */
	typedef dasp::Array<Assignment> FrameAssignment;

	/** An edge connecting to clusters */
	struct Edge
	{
		cluster_id_type a, b;

		friend bool operator<(const Edge& x, const Edge& y) {
			return x.a < y.a || (x.a == y.a && x.b < y.b);
		}
	};

	/** Various data related to one timestamp */
	struct Frame
	{
		int time;
		RgbdData rgbd;
		ClusterList clusters;
		FrameAssignment assignment;
		std::vector<Edge> edges;
	};

	typedef std::shared_ptr<Frame> FramePtr;

	/** Creates a frame from color and depht data */
	RgbdData CreateRgbdData(const Rgbd& data);

	/** Computes cluster density */
	Eigen::MatrixXf ComputeClusterDensity(const RgbdData& rgbd);

	/** Samples clusters from cluster density */
	std::vector<Cluster> SampleClustersFromDensity(const RgbdData& rgbd, const Eigen::MatrixXf& density);

	/** Creates a frame from RGBD data and initial clusters */
	FramePtr CreateFrame(int time, const RgbdData& rgbd, const std::vector<Cluster>& clusters);

	/** A timeseries of frames */
	struct Timeseries
	{
		std::vector<FramePtr> frames;

		ClusterContainer clusters;

		int slices() const {
			return frames.size();
		}

		int rows() const {
			return frames.front()->rgbd.rows();
		}

		int cols() const {
			return frames.front()->rgbd.cols();
		}

		int getBeginTime() const {
			return frames.empty() ? 0 : frames.front()->time;
		}

		int getEndTime() const {
			return frames.empty() ? 0 : frames.back()->time+1;
		}

		int getDuration() const {
			return getEndTime() - getBeginTime();
		}

		const FramePtr& getFrame(int t) const {
			const int t0 = getBeginTime();
			assert(t0 <= t && t < t0 + frames.size());
			return frames[t-t0];
		}

		/** Get all frames with ta <= frame_time <= tb */
		std::vector<FramePtr> getFrameRange(int ta, int tb) const {
			assert(ta <= tb);
			int i1 = 0;
			while(i1 < frames.size() && frames[i1]->time < ta) i1++;
			int i2 = i1;
			while(i2 < frames.size() && frames[i2]->time <= tb) i2++;
//			std::cout << t_begin << " " << t_end << " -> " << i1 << " " << i2 << std::endl;
			return std::vector<FramePtr>(frames.begin() + i1, frames.begin() + i2);
		}

		void add(const FramePtr& f) {
			frames.push_back(f);
		}

		std::vector<FramePtr> purge(int tmin) {
			int i;
			for(i=0; i<frames.size() && frames[i]->time < tmin; ++i);
			std::vector<FramePtr> purged_frames(frames.begin(), frames.begin() + i);
			frames.erase(frames.begin(), frames.begin() + i);
			return purged_frames;
		}

	};

	/** Updates clusters near a given timestep */
	void UpdateClusters(int time, Timeseries& timeseries);

	/** Graph structure for weighted cluster graph */
	typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS,
		Cluster, // vertex annotation
		boost::property<boost::edge_weight_t, float> // edge annotation
	> ClusterGraph;

	/** Computes cluster graph of frames */
	ClusterGraph ComputeClusterGraph(const std::vector<FramePtr>& frames);

	/** Writes clusters to a file */
	void IOWriteClusters(const std::string& fn, const std::vector<Cluster>& clusters);

	/** Reads clusters from a file */
	std::vector<Cluster> IOReadClusters(const std::string& fn);

	/** Writes edges to a file */
	void IOWriteEdges(const std::string& fn, const std::vector<Edge>& edges);

	/** Writes a cluster graph to files */ 
	void IOWriteGraph(const std::string& fn_vertices, const std::string& fn_edges, const ClusterGraph& graph);

	/** Reads a cluster graph from files */ 
	ClusterGraph IOReadGraph(const std::string& fn_vertices, const std::string& fn_edges);

	/** Control structure for continuous temporal-depth-adative superpoint generation */
	struct ContinuousSupervoxels
	{
		/** Starts temporal superpoint generation */
		void start();

		/** Processes on timestamp and updates/generates temporal superpoint */
		void step(const Rgbd& rgbd);

		/** Gets number of active clusters */
		int numActiveClusters() const;

		/** Gets number of inactive clusters */
		int numInactiveClusters() const;

		/** Gets total number of clusters */
		int numClusters() const;

		/** Gets all clusters */
		std::vector<Cluster> getAllClusters() const;

		/** Gets total cluster graph */
		const ClusterGraph& getGraph() const {
			return graph_;
		}

	private:
		bool is_first_;
		Eigen::MatrixXf last_density_;
		Timeseries series_;
		ClusterGraph graph_;
	};

	enum class PlotStyle {
		Color, Age, AssignmentDistance, Label, ClusterBorder
	};

	/** Plots superpixel cluster color */
	slimage::Image3ub DebugPlotClusters(const FramePtr& frame, const std::vector<PlotStyle>& styles);

	/** Computes compression error
	 *   sum_i (mu_{a(i)} - mu)^2 / sum_i (x_i - mu)^2
	 * where
	 *   x_i value of i-th pixel,
	 *   mu_j mean for j-th cluster,
	 *   a(i) cluster assignment of i-th pixel,
	 *   mu total mean of all pixels
	 */
	Eigen::Vector2f EvaluateComputeCompressionError(const FramePtr& frame);

	/** Computes the compression error for a downsampling with approx the same number of clusters */
	Eigen::Vector2f EvaluateComputeDownsampleCompressionError(const FramePtr& frame);

}

#endif
