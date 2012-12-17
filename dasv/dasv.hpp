#ifndef DASV_DASV_HPP
#define DASV_DASV_HPP

#include <Slimage/Slimage.hpp>
#include <Eigen/Dense>

namespace dasv
{
	struct Point
	{
		Eigen::Vector3f color;
		Eigen::Vector3f position;
		Eigen::Vector3f normal;
		float cluster_radius_px;
		bool valid;
	};

	struct Frame
	{
		typedef std::vector<Point> Container;
		typedef Container::iterator it;
		typedef Container::const_iterator cit;

		int rows, cols;
		Container points;

		Frame()
		: rows(0), cols(0)
		{}

		Frame(int nrows, int ncols)
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
			return (0 <= i && i < rows && 0 <= j && j < cols);
		}

		const Point& operator()(int i, int j) const {
			return points[i*cols+j];
		}

		Point& operator()(int i, int j) {
			return points[i*cols+j];
		}
		
	};

	Frame CreateFrame(const slimage::Image3ub& color, const slimage::Image1ui16& depth);

	Eigen::MatrixXf ComputeClusterDensity(const Frame& frame);

	struct Cluster
	{
		Eigen::Vector3f color;
		Eigen::Vector3f position;
		Eigen::Vector3f normal;
		float time;
	};

	std::vector<Cluster> SampleClusters(const Eigen::MatrixXf& density, float time);

}

#endif
