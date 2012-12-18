#ifndef DASV_DASV_HPP
#define DASV_DASV_HPP

#include <Slimage/Slimage.hpp>
#include <Eigen/Dense>
#include <vector>
#include <map>
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
	};

	/** Samples clusters from cluster density */
	std::vector<Cluster> SampleClustersFromDensity(const RgbdData& frame, const Eigen::MatrixXf& density);

	/** Various data related to one timestamp */
	struct Frame
	{
		RgbdData rgbd;
		std::vector<Cluster> clusters;
		int time;
	};

	typedef std::shared_ptr<Frame> FramePtr;

	/** A series of frames and clusters */
	struct Timeseries
	{
		std::vector<FramePtr> frames;
		int t0;

		const FramePtr& getFrame(int t) {
			assert(t0 <= t && t < t0 + frames.size());
			return frames[t-t0];
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
