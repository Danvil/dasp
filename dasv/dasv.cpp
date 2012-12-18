/**
 * Code conventions:
 * - pixel indices are (signed) INT
 * - linear algebra is used via EIGEN
 * - Eigen matrices storage order is COLUMN-MAJOR
 *   The following loop structure should be used:
 *   for(int j=0; j<cols; j++) {
 *     for(int i=0; i<rows; i++) {
 *	     m(i,j) = 1.0f;
 *   }}
 * - With (x,y) coordinates these correspondences should be used:
 *     rows -> width
 *     cols -> height
 *     m(x,y) (it is optimal to use x in the inner loop)
 */

#include "dasv.hpp"
#include <random>
#include <assert.h>

namespace dasv
{

std::mt19937 random_engine;

constexpr float DEPTH_TO_Z = 0.001f;
constexpr float CENTER_X = 320.0f;
constexpr float CENTER_Y = 240.0f;
constexpr float PX_FOCAL = 528.0f;
constexpr float CLUSTER_RADIUS = 0.03f;
constexpr float CLUSTER_TIME = 30.0f; // 1 second

constexpr float PI = 3.1415f;

void ComputeRgbdDataNormals(RgbdData& frame)
{
	const int rows = frame.rows;
	const int cols = frame.cols;
	for(int y=0; y<cols; y++) {
		for(int x=0; x<rows; x++) {
			frame(x,y).normal = Eigen::Vector3f(0,0,-1); // FIXME implement
		}
	}
}

RgbdData CreateRgbdData(const slimage::Image3ub& img_color, const slimage::Image1ui16& img_depth)
{
	const int rows = img_color.height();
	const int cols = img_color.width();
	RgbdData f(rows, cols);
	for(int y=0, i=0; y<cols; y++) {
		for(int x=0; x<rows; x++, i++) {
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
			// normal -> ComputeRgbdDataNormals
			// world cluster radius
			point.cluster_radius_px = CLUSTER_RADIUS / z_over_f;
			// valid
			point.valid = (depth != 0);
		}
	}
	ComputeRgbdDataNormals(f);
	return f;
}

Eigen::MatrixXf ComputeClusterDensity(const RgbdData& frame)
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

bool FindValidSeedPoint(const RgbdData& points, int& sx0, int& sy0, int range)
{
	if(range == 0) {
		return 0 <= sx0 && sx0 < points.rows
			&& 0 <= sy0 && sy0 < points.cols
			&& points(sx0,sy0).valid;
	}
	// add random offset to add noise
//	std::uniform_int<int> delta(-range, +range); // FIXME
	unsigned int trials = 0;
	while(trials < 100) {
		int sx = sx0;// + delta(random_engine);
		int sy = sy0;// + delta(random_engine);
		if(    0 <= sx && sx < points.rows
			&& 0 <= sy && sy < points.cols
			&& points(sx,sy).valid
		) {
			sx0 = sx;
			sy0 = sy;
			return true;
		}
		trials++;
	}
	return false;
}

void SampleDensityImplRec(
		const RgbdData& points,
		std::vector<Eigen::Vector2f>& seeds,
		const std::vector<Eigen::MatrixXf>& mipmaps,
		std::vector<Eigen::MatrixXf>& carry_mipmaps,
		int level, int x, int y)
{
	const Eigen::MatrixXf& mm = mipmaps[level];
	Eigen::MatrixXf& carry_mm = carry_mipmaps[level];

	// compute density by multiplying percentage with parent total
	float v = mm(x, y) + carry_mm(x, y);

	// FIXME if low density is carried over on a low-res mipmap
	// FIXME and the target cell has a high density
	// FIXME the carried over density is not considered on a high-res mipmap level

	if(level <= 1 || v <= 1.5f) {
		if(v >= 0.5f) {
			// set seed point in the middel of the cell
			int half = (1 << (level - 1));
			int sx = (x << level) + half;
			int sy = (y << level) + half;
			if(FindValidSeedPoint(points, sx, sy, half/4)) { // place near center
				seeds.push_back(Eigen::Vector2f(sx, sy));
				// reduce density by 1
				v -= 1.0f;
			}
		}
		// distribute remaining density to neighbours
		// mm(x+1,y  ) += 7.0f / 16.0f * v;
		// mm(x-1,y+1) += 3.0f / 16.0f * v;
		// mm(x  ,y+1) += 5.0f / 16.0f * v;
		// mm(x+1,y+1) += 1.0f / 16.0f * v;
		// with range test *sigh*
		float q = 0.0f;
		bool xm1ok = (0 < x);
		bool xp1ok = (x+1 < mm.rows());
		bool yp1ok = (y+1 < mm.cols());
		if(xp1ok) 			q += 7.0f;
		if(yp1ok) {
			if(xm1ok) 		q += 3.0f;			
							q += 5.0f;
			if(xp1ok) 		q += 1.0f;
		}
		if(q > 0) {
			float scl = v / q;
			if(xp1ok) 		carry_mm(x+1,y  ) += 7.0f * scl;
			if(yp1ok) {
				if(xm1ok) 	carry_mm(x-1,y+1) += 3.0f * scl;			
							carry_mm(x  ,y+1) += 5.0f * scl;
				if(xp1ok) 	carry_mm(x+1,y+1) += 1.0f * scl;
			}
		}
	}
	else {
		// go down
		SampleDensityImplRec(points, seeds, mipmaps, carry_mipmaps, level - 1, 2*x,     2*y    );
		SampleDensityImplRec(points, seeds, mipmaps, carry_mipmaps, level - 1, 2*x,     2*y + 1);
		SampleDensityImplRec(points, seeds, mipmaps, carry_mipmaps, level - 1, 2*x + 1, 2*y    );
		SampleDensityImplRec(points, seeds, mipmaps, carry_mipmaps, level - 1, 2*x + 1, 2*y + 1);
	}
}

int find_next_pow2(int x)
{
	int a = 1;
	while(x < a) a *= 2;
	return a;
}

Eigen::MatrixXf ComputeMipmap(const Eigen::MatrixXf& data)
{
	const int rows = data.rows();
	const int cols = data.cols();
	assert(rows % 2 == 0);
	assert(cols % 2 == 0);
	const int mm_rows = find_next_pow2(rows)/2;
	const int mm_cols = find_next_pow2(cols)/2;
	Eigen::MatrixXf mm = Eigen::MatrixXf::Zero(mm_rows, mm_cols);
	for(int y=0; y<cols; y+=2) {
		for(int x=0; x<rows; x+=2) {
			float q = data(x,y) + data(x,y+1) + data(x+1,y) + data(x+1,y+1);
			mm(x/2, y/2) = q;
		}
	}
	return mm;
}

std::vector<Eigen::MatrixXf> ComputeMipmaps(const Eigen::MatrixXf& data, int min_size)
{
	std::vector<Eigen::MatrixXf> mm;
	mm.reserve(10); // 2^10 = 1024
	mm.push_back(data);
	while(true) {
		const Eigen::MatrixXf& q = mm.back();
		if(q.rows() <= min_size || q.cols() <= min_size) {
			break;
		}
		mm.push_back(ComputeMipmap(q));
	}
	return mm;
}

std::vector<Eigen::Vector2f> SampleDensityImpl(const RgbdData& frame, const Eigen::MatrixXf& density)
{
	// compute mipmaps
	std::vector<Eigen::MatrixXf> mipmaps = ComputeMipmaps(density, 32);
	std::vector<Eigen::MatrixXf> carry_mipmaps(mipmaps.size());
	for(unsigned int i=1; i<mipmaps.size(); i++) { // HACK: skip first as we never use it
		carry_mipmaps[i] = Eigen::MatrixXf::Zero(mipmaps[i].rows(), mipmaps[i].cols());
	}
	// now create pixel seeds
	std::vector<Eigen::Vector2f> seeds;
	seeds.reserve(1000);
	SampleDensityImplRec(frame, seeds, mipmaps, carry_mipmaps, mipmaps.size() - 1, 0, 0);
	return seeds;
}

std::vector<Cluster> SampleClustersFromDensity(const RgbdData& frame, const Eigen::MatrixXf& density)
{
	std::vector<Eigen::Vector2f> points = SampleDensityImpl(frame, density);
	// clusters 
	std::vector<Cluster> clusters(points.size());
	for(int i=0; i<clusters.size(); i++) {
		Eigen::Vector2f px = points[i];
		int x = static_cast<int>(px.x());
		int y = static_cast<int>(px.y());
		if(!frame.valid(x,y)) {
			// skip
			continue;
		}
		const Point& fp = frame(y,x);
		Cluster& c = clusters[i];
		c.pixel = px;
		// FIXME use small neighbourhood
		c.color = fp.color;
		c.position = fp.position;
		c.normal = fp.normal;
	}
	return clusters;
}

void UpdateClusters(int time, Timeseries& timeseries)
{

}

}
