#ifndef DASV_DASV_HPP
#define DASV_DASV_HPP

#include <Slimage/Slimage.hpp>
#include <Eigen/Dense>
#include <boost/multi_array.hpp>
#include <vector>
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

	/** Various data related to one timestamp */
	struct Frame
	{
		int time;
		RgbdData rgbd;
		std::vector<ClusterPtr> clusters;
		FrameAssignment assignment;
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
			result.insert(result.begin(), frames.begin() + i1, frames.begin() + i2);
			return result;
		}

		void add(const FramePtr& f) {
			frames.push_back(f);
		}

		std::vector<ClusterPtr> purge(int tmin) {
			std::vector<ClusterPtr> clusters;
			int i;
			for(i=0; i<frames.size() && frames[i]->time < tmin; ++i) {
				clusters.insert(clusters.begin(), frames[i]->clusters.begin(), frames[i]->clusters.end());
			}
			frames.erase(frames.begin(), frames.begin() + i);
			return clusters;
		}

	};

	/** Updates clusters near a given timestep */
	void UpdateClusters(int time, Timeseries& timeseries);

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

		/** Gets all clusters */
		std::vector<Cluster> getAllClusters() const;

	private:
		bool is_first_;
		Eigen::MatrixXf last_density_;
		Timeseries series_;
		std::vector<ClusterPtr> inactive_clusters_;
	};

	/** Opens a window to display matrix data */
	void DebugShowMatrix(const std::string& filename, const Eigen::MatrixXf& mat, float scl);

	/** Writes clusters to a file */
	void DebugWriteClusters(const std::string& fn, const std::vector<Cluster>& clusters);

	/** Computes superpixel image from clusters and assignment */
	slimage::Image3ub DebugCreateSuperpixelImage(const FramePtr& frame);

	/** Computes compression error
	 *   sum_i (mu_{a(i)} - mu)^2 / sum_i (x_i - mu)^2
	 * where
	 *   x_i value of i-th pixel,
	 *   mu_j mean for j-th cluster,
	 *   a(i) cluster assignment of i-th pixel,
	 *   mu total mean of all clusters
	 */
	Eigen::Vector2f EvaluateComputeCompressionError(const FramePtr& frame);

}

#endif
