#ifndef DASV_DASV_HPP
#define DASV_DASV_HPP

#include <Slimage/Slimage.hpp>
#include <Eigen/Dense>
#include <boost/graph/adjacency_list.hpp>
#include <vector>
#include <ostream>
#include <string>
#include <memory>

namespace dasv
{
	/** 2D data structure implemented with an std::vector */
	template<typename T>
	struct Vector2D
	{
	public:
		typedef std::vector<T> Container;
		typedef typename Container::iterator it;
		typedef typename Container::const_iterator const_it;
		Vector2D() : rows_(0), cols_(0) {}
		Vector2D(int rows, int cols) : rows_(rows), cols_(cols), data_(rows*cols) { }
		Vector2D(int rows, int cols, const T& t) : rows_(rows), cols_(cols), data_(rows*cols, t) { }
		Container& data() { return data_; }
		const Container& data() const { return data_; }
		it begin() { return data_.begin(); }
		const_it begin() const { return data_.begin(); }
		it end() { return data_.end(); }
		const_it end() const { return data_.end(); }
		int size() const { return rows_*cols_; }
		T& operator[](int i) { return data_[i]; }
		const T& operator[](int i) const { return data_[i]; }
		int rows() const { return rows_; }
		int cols() const { return cols_; }
		bool isValid(int i, int j) const { return (0 <= i && i < rows_ && 0 <= j && j < cols_); }
		T& operator()(int i, int j) { return data_[i + j*rows_]; }
		const T& operator()(int i, int j) const { return data_[i + j*rows_]; }
	private:
		int rows_, cols_;
		std::vector<T> data_;
	};

	/** Voxel data */
	struct Point
	{
		Eigen::Vector3f color;
		Eigen::Vector3f position;
		Eigen::Vector3f normal;
		float cluster_radius_px;
		bool valid;
	};

	/** A frame of voxels at a given timestamp */
	typedef Vector2D<Point> RgbdData;

	/** A voxel cluster (supervoxel) */
	struct Cluster
	{
		Eigen::Vector2f pixel;
		Eigen::Vector3f color;
		Eigen::Vector3f position;
		Eigen::Vector3f normal;
		float cluster_radius_px;
		int time;
		int id;
		bool valid;

		int unique_id() const {
			return time * (1<<16) + id;
		}

		friend std::ostream& operator<<(std::ostream& os, const Cluster& c)
		{
			os << c.time
				<< "\t" << c.id
				<< "\t" << (c.valid ? 1 : 0)
				<< "\t" << c.cluster_radius_px
				<< "\t" << c.pixel.x() << "\t" << c.pixel.y()
				<< "\t" << c.color.x() << "\t" << c.color.y() << "\t" << c.color.z()
				<< "\t" << c.position.x() << "\t" << c.position.y() << "\t" << c.position.z()
				<< "\t" << c.normal.x() << "\t" << c.normal.y() << "\t" << c.normal.z();
			return os;
		}

	};

	typedef std::shared_ptr<Cluster> ClusterPtr;

	/** Pixel to cluster assignment */
	struct Assignment
	{
		ClusterPtr cluster;
		float distance;
	};

	/** Frame pixel to cluster assignment */
	typedef Vector2D<Assignment> FrameAssignment;

	/** An edge connecting to clusters */
	struct Edge
	{
		ClusterPtr a, b;

		int id() const {
			return a->unique_id() | b->unique_id();
		}
		
		friend bool operator==(const Edge& x, const Edge& y) {
			return (x.a == y.a && x.b == y.b) || (x.a == y.b && x.b == y.a);
		}

		friend bool operator<(const Edge& x, const Edge& y) {
			return x.id() < y.id();
		}
	};

	/** Various data related to one timestamp */
	struct Frame
	{
		int time;
		RgbdData rgbd;
		std::vector<ClusterPtr> clusters;
		FrameAssignment assignment;
		std::vector<Edge> edges;
	};

	typedef std::shared_ptr<Frame> FramePtr;

	/** Creates a frame from color and depht data */
	RgbdData CreateRgbdData(const slimage::Image3ub& color, const slimage::Image1ui16& depth);

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

		std::vector<FramePtr> getFrameRange(int t_begin, int t_end) const {
			assert(t_begin < t_end);
			int i1 = 0;
			while(i1 < frames.size() && frames[i1]->time < t_begin) i1++;
			int i2 = i1;
			while(i2 < frames.size() && frames[i2]->time < t_end) i2++;
			std::vector<FramePtr> result;
//			std::cout << t_begin << " " << t_end << " -> " << i1 << " " << i2 << std::endl;
			result.insert(result.begin(), frames.begin() + i1, frames.begin() + i2);
			return result;
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

	namespace detail
	{
		template<class Iter> struct iter_pair_range : std::pair<Iter,Iter>
		{
			iter_pair_range(const std::pair<Iter,Iter>& x) : std::pair<Iter,Iter>(x) {}
			Iter begin() const { return this->first; }
			Iter end() const { return this->second; }
		};
	}

	template<class Iter>
	inline detail::iter_pair_range<Iter> as_range(const std::pair<Iter,Iter>& x) {
		return detail::iter_pair_range<Iter>(x);
	}

	/** Computes cluster graph of frames */
	ClusterGraph ComputeClusterGraph(const std::vector<FramePtr>& frames);

	/** Writes clusters to a file */
	void IOWriteClusters(const std::string& fn, const std::vector<ClusterPtr>& clusters);

	/** Writes edges to a file */
	void IOWriteEdges(const std::string& fn, const std::vector<Edge>& edges);

	/** Writes a cluster graph to files */ 
	void IOWriteGraph(const std::string& fn_vertices, const std::string& fn_edges, const ClusterGraph& graph);

	/** Control structure for continuous temporal-depth-adative superpoint generation */
	struct ContinuousSupervoxels
	{
		/** Starts temporal superpoint generation */
		void start(int rows, int cols);

		/** Processes on timestamp and updates/generates temporal superpoint */
		void step(const slimage::Image3ub& color, const slimage::Image1ui16& depth);

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
		std::vector<ClusterPtr> inactive_clusters_;
		ClusterGraph graph_;
		std::vector<Edge> delayed_edges_;
	};

	/** Opens a window to display matrix data */
	void DebugShowMatrix(const std::string& filename, const Eigen::MatrixXf& mat, float scl);

	/** Computes superpixel image from clusters and assignment */
	slimage::Image3ub DebugCreateSuperpixelImage(const FramePtr& frame, bool borders, bool age_colors);

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
