#ifndef DASV_DASV_HPP
#define DASV_DASV_HPP

#include <Slimage/Slimage.hpp>
#include <Eigen/Dense>
#include <boost/multi_array.hpp>
#include <vector>
#include <memory>

namespace dasv
{
	/** Voxel data */
	struct Point
	{
		Eigen::Vector3f color;
		Eigen::Vector3f position;
		Eigen::Vector3f normal;
		float cluster_radius_px;
		bool valid;
	};

	inline bool IsValidIndex(int rows, int cols, int i, int j) {
		return (0 <= i && i < rows && 0 <= j && j < cols);
	}

	/** A frame of voxels at a given timestamp
	 * COLUMN-MAJOR storage order!
	 */
	struct RgbdData
	{
		typedef std::vector<Point> Container;
		typedef Container::iterator it;
		typedef Container::const_iterator cit;

		int rows, cols;
		Container points;

		RgbdData()
		: rows(0), cols(0)
		{}

		RgbdData(int nrows, int ncols)
		: rows(nrows), cols(ncols),
		  points(nrows*ncols)
		{}

		int size() const { return rows*cols; }

		cit begin() const { return points.begin(); }

		it begin() { return points.begin(); }

		cit end() const { return points.end(); }

		it end() { return points.end(); }

		const Point& operator[](int i) const {
			return points[i];
		}

		Point& operator[](int i) {
			return points[i];
		}

		bool valid(int i, int j) const {
			return IsValidIndex(rows, cols, i, j);
		}

		const Point& operator()(int i, int j) const {
			return points[i+j*rows];
		}

		Point& operator()(int i, int j) {
			return points[i+j*rows];
		}
		
	};

	/** Creates a frame from color and depht data */
	RgbdData CreateRgbdData(const slimage::Image3ub& color, const slimage::Image1ui16& depth);

	/** Computes cluster density */
	Eigen::MatrixXf ComputeClusterDensity(const RgbdData& frame);

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

	/** Samples clusters from cluster density */
	std::vector<Cluster> SampleClustersFromDensity(const RgbdData& frame, const Eigen::MatrixXf& density);

	struct Assignment
	{
		ClusterPtr cluster;
		float distance;
	};

	/** Frame pixel to cluster assignment */
	typedef boost::multi_array<Assignment, 2> assigment_type;

	/** Various data related to one timestamp */
	struct Frame
	{
		int time;
		RgbdData rgbd;
		std::vector<ClusterPtr> clusters;
		assigment_type assignment;
	};

	typedef std::shared_ptr<Frame> FramePtr;

	/** Creates a frame from RGBD data and initial clusters */
	FramePtr CreateFrame(int time, const RgbdData& rgbd, const std::vector<Cluster>& clusters);

	/** A series of frames and clusters */
	struct Timeseries
	{
		std::vector<FramePtr> frames;
		int t0;

		const FramePtr& getFrame(int t) {
			assert(t0 <= t && t < t0 + frames.size());
			return frames[t-t0];
		}

		std::vector<FramePtr> getFrameRange(int t1, int t2) const {
			assert(t0 <= t1 && t2 < t0+frames.size());
			std::vector<FramePtr> result;
			result.insert(result.begin(), frames.begin() + t1 - t0, frames.begin() + t2 - t0);
			return result;
		}
	};

	/**
	 * Updates clusters for a given timestep
	 * @param time current timestep
	 * @param frames time -> frames
	 * @param clusters time -> clusters for this timestep
	 */
	void UpdateClusters(int time, Timeseries& timeseries);

}

#endif
