#include "dasv.hpp"

namespace dasv
{

constexpr float DEPTH_TO_Z = 0.001f;
constexpr float CENTER_X = 320.0f;
constexpr float CENTER_Y = 240.0f;
constexpr float PX_FOCAL = 528.0f;
constexpr float CLUSTER_RADIUS = 0.03f;
constexpr float CLUSTER_TIME = 30.0f; // 1 second

constexpr float PI = 3.1415f;

void ComputeFrameNormals(Frame& frame)
{
	const int rows = frame.rows;
	const int cols = frame.cols;
	for(int y=0; y<rows; y++) {
		for(int x=0; x<cols; x++) {
			frame(y,x).normal = Eigen::Vector3f(0,0,-1); // FIXME implement
		}
	}
}

Frame CreateFrame(const slimage::Image3ub& img_color, const slimage::Image1ui16& img_depth)
{
	const int rows = img_color.height();
	const int cols = img_color.width();
	Frame f(rows, cols);
	for(int y=0, i=0; y<rows; y++) {
		for(int x=0; x<cols; x++, i++) {
			Point& point = f[i];
			const uint16_t depth = img_depth[i];
			const float z_over_f = DEPTH_TO_Z * static_cast<float>(depth) / PX_FOCAL;
			// RGB color
			const slimage::Pixel3ub& color = img_color[i];
			point.color = (1.0f/255.0f) * Eigen::Vector3f(
					static_cast<float>(color[0]),
					static_cast<float>(color[1]),
					static_cast<float>(color[2]));
			// point from depth
			point.position = z_over_f * Eigen::Vector3f(
					static_cast<float>(x) - CENTER_X,
					static_cast<float>(y) - CENTER_Y,
					PX_FOCAL);
			// normal -> ComputeFrameNormals
			// world cluster radius
			point.cluster_radius_px = CLUSTER_RADIUS / z_over_f;
			// valid
			point.valid = (depth != 0);
		}
	}
	ComputeFrameNormals(frame);
	return f;
}

Eigen::MatrixXf ComputeClusterDensity(const Frame& frame)
{
	Eigen::MatrixXf density(frame.rows, frame.cols);
	const int size = frame.size();
	for(int i=0; i<size; i++) {
		const Point& p = frame[i];
		// rho = r_px|^2 * pi / sqrt(||g||^2+1)
		// 1/sqrt(||g||^2+1) = n_z because g = -(n_x/n_z, n_y/n_z)
		// TODO n_z should always be negative so abs(n_z) should equal -n_z.
		const float A = p.cluster_radius_px * p.cluster_radius_px * PI * std::abs(p.normal.z());
		const float rho = 1.0f / A / CLUSTER_TIME;
		density(i) = rho; // FIXME is density(i) correct?
	}
	return density;
}


std::vector<Eigen::Vector2f> Sample(const Eigen::MatrixXf& density)
{
	std::vector<Eigen::Vector2f> points;
	points.reserve(1000);
	return points;
}

std::vector<Cluster> SampleClusters(const Eigen::MatrixXf& density, float time)
{
	std::vector<Eigen::Vector2f> points = Sample(density);
	// clusters 
	// std::vector<Cluster> clusters()
	// cluster.reserve(1000);

	return clusters;
}


}
