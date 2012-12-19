/**
 * Code conventions:
 * - pixel indices are (signed) INT
 * - Floating point matrices and linear algebra use Eigen
 *   Eigen matrices storage order is COLUMN-MAJOR
 *   The following loop structure should be used:
 *     Eigen::MatrixXf m(rows, cols);
 *     for(int i=0; i<cols; i++)
 *       for(int j=0; j<rows; j++)
 *	       m(j,i) = 42.0f;
 * - User defined 2D arrays use boost::multi_array
 *   boost::multi_array storage order is C-STYLE
 *   The following loop structure should be used:
 *     boost::multi_array<T> a(boost::extents[cols][rows]);
 *     for(int i=0; i<a.shape()[0]; i++)
 *       for(int j=0; j<a.shape()[1]; j++)
 *	       a[i][j] = 42.0f;
 * - With (x,y) coordinates these correspondences should be used:
 *     width -> rows -> a.shape()[1]
 *     height -> cols -> a.shape()[0]
 *     m(x,y) (it is optimal to use x in the inner loop)
 *     a[y][x] (it is optimal to use x in the inner loop)
 */

#include "dasv.hpp"
#include <Slimage/Gui.hpp>
#include <boost/multi_array.hpp>
#include <boost/lexical_cast.hpp>
#include <random>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <assert.h>

//#define GUI_DEBUG_VERBOSE
//#define GUI_DEBUG_NORMAL

namespace dasv
{

std::mt19937 random_engine;

constexpr float DEPTH_TO_Z = 0.001f;
constexpr float CENTER_X = 320.0f;
constexpr float CENTER_Y = 240.0f;
constexpr float PX_FOCAL = 528.0f;
constexpr float CLUSTER_RADIUS = 0.025f;
constexpr int CLUSTER_TIME_RADIUS = 5; // TR=15 -> 0.5 s
constexpr int CLUSTER_ITERATIONS = 1;
constexpr float CLUSTER_RADIUS_MULT = 1.5f;
constexpr float SPATIAL_TIME_INCREASE = 0.005f; // 0.15 m/s

constexpr float PI = 3.1415f;

void DebugShowMatrix(const std::string& filename, const Eigen::MatrixXf& mat, float scl)
{
#ifdef GUI_DEBUG_NORMAL
	slimage::Image1f img(mat.rows(), mat.cols());
	img.buffer().copyFrom(&mat(0,0));
	slimage::gui::Show(filename, img, scl, 0);
#endif
}

Eigen::MatrixXf DebugDoubleMatrixSize(const Eigen::MatrixXf& mat, int n)
{
	Eigen::MatrixXf last = mat;
	Eigen::MatrixXf result;
	for(int k=0; k<n; k++) {
		result = Eigen::MatrixXf(last.rows()*2, last.cols()*2);
		for(int i=0; i<result.cols(); i++)
			for(int j=0; j<result.rows(); j++)
				result(j,i) = last(j/2,i/2);
		last = result;
	}
	return last;
}

void ComputeRgbdDataNormals(RgbdData& rgbd)
{
	const int NY = rgbd.cols();
	const int NX = rgbd.rows();
	for(int y=0; y<NY; y++) {
		for(int x=0; x<NX; x++) {
			rgbd(x,y).normal = Eigen::Vector3f(0,0,-1); // FIXME implement
		}
	}
}

Eigen::Vector2f CameraProject(const Eigen::Vector3f& p)
{
	return {
		CENTER_X + PX_FOCAL*p.x()/p.z(),
		CENTER_Y + PX_FOCAL*p.y()/p.z()
	};
}

RgbdData CreateRgbdData(const slimage::Image3ub& img_color, const slimage::Image1ui16& img_depth)
{
	const int NX = img_color.width();
	const int NY = img_color.height();
	RgbdData rgbd(NX, NY);
	for(int y=0, i=0; y<NY; y++) {
		for(int x=0; x<NX; x++, i++) {
			Point& point = rgbd(x,y);
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
	ComputeRgbdDataNormals(rgbd);
	return rgbd;
}

Eigen::MatrixXf ComputeFrameDensity(const RgbdData& rgbd)
{
	const int NY = rgbd.cols();
	const int NX = rgbd.rows();
	Eigen::MatrixXf density(NX, NY);
	for(int y=0; y<NY; y++) {
		for(int x=0; x<NX; x++) {
			const Point& p = rgbd(x,y);
			if(p.valid) {
				// rho = r_px|^2 * pi / sqrt(||g||^2+1)
				// 1/sqrt(||g||^2+1) = n_z because g = -(n_x/n_z, n_y/n_z)
				// TODO n_z should always be negative so abs(n_z) should equal -n_z.
				const float A = p.cluster_radius_px * p.cluster_radius_px * PI * std::abs(p.normal.z());
				const float rho = 1.0f / A;
				density(x,y) = rho; // FIXME is density(i) correct?
			}
			else {
				density(x,y) = 0.0f;
			}
		}
	}
	return density;
}

Eigen::MatrixXf ComputeClusterDensity(int rows, int cols, const std::vector<Cluster>& clusters)
{
	Eigen::MatrixXf density = Eigen::MatrixXf::Zero(rows, cols);
	// range R of kernel is s.t. phi(x) >= 0.01 * phi(0) for all x <= R
	const float cRange = 1.21f; // BlueNoise::KernelFunctorInverse(0.01f);
	for(const Cluster& c : clusters) {
		if(!c.valid) {
			continue;
		}
		const int sx = static_cast<int>(c.pixel.x() + 0.5f);
		const int sy = static_cast<int>(c.pixel.y() + 0.5f);
		const float sxf = c.pixel.x();
		const float syf = c.pixel.y();
		// seed corresponds to a kernel at position (x,y) with sigma = rho(x,y)^(-1/2)
		const float rho = 1.0f / (c.cluster_radius_px*c.cluster_radius_px*PI*std::abs(c.normal.z()));
		// kernel influence range
		const int R = static_cast<int>(std::ceil(cRange / std::sqrt(rho)));
		const int xmin = std::max<int>(sx - R, 0);
		const int xmax = std::min<int>(sx + R, rows - 1);
		const int ymin = std::max<int>(sy - R, 0);
		const int ymax = std::min<int>(sy + R, cols - 1);
		for(int yi=ymin; yi<=ymax; yi++) {
			for(int xi=xmin; xi<=xmax; xi++) {
				const float dx = static_cast<float>(xi) - sxf;
				const float dy = static_cast<float>(yi) - syf;
				const float d2 = dx*dx + dy*dy;
				const float delta = rho * std::exp(-PI*rho*d2);// BlueNoise::KernelFunctorSquare(rho*d2);
				density(xi, yi) += delta;
			}
		}
	}
	return density;
}

bool FindValidSeedPoint(const RgbdData& rgbd, int& sx0, int& sy0, int range)
{
	const int NY = rgbd.cols();
	const int NX = rgbd.rows();
	if(range == 0) {
		return 0 <= sx0 && sx0 < NX
			&& 0 <= sy0 && sy0 < NY
			&& rgbd(sx0,sy0).valid;
	}
	// add random offset to add noise
	std::uniform_int_distribution<int> delta(-range, +range); // FIXME
	unsigned int trials = 0;
	while(trials < 100) {
		int sx = sx0 + delta(random_engine);
		int sy = sy0 + delta(random_engine);
		if(    0 <= sx && sx < NX
			&& 0 <= sy && sy < NY
			&& rgbd(sx,sy).valid
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
		const RgbdData& rgbd,
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
			if(FindValidSeedPoint(rgbd, sx, sy, half/2)) { // place near center
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
		SampleDensityImplRec(rgbd, seeds, mipmaps, carry_mipmaps, level - 1, 2*x,     2*y    );
		SampleDensityImplRec(rgbd, seeds, mipmaps, carry_mipmaps, level - 1, 2*x,     2*y + 1);
		SampleDensityImplRec(rgbd, seeds, mipmaps, carry_mipmaps, level - 1, 2*x + 1, 2*y    );
		SampleDensityImplRec(rgbd, seeds, mipmaps, carry_mipmaps, level - 1, 2*x + 1, 2*y + 1);
	}
}

inline int find_next_pow2(int x)
{
	int a = 1;
	while(a < x) a *= 2;
	return a;
}

Eigen::MatrixXf ComputeMipmap(const Eigen::MatrixXf& data)
{
	const int rows = data.rows();
	const int cols = data.cols();
//	std::cout << "high: " << rows << " " << cols << std::endl;
	assert(rows % 2 == 0);
	assert(cols % 2 == 0);
	const int size = find_next_pow2(std::max(rows,cols))/2;
//	std::cout << "low: " << mm_rows << " " << mm_cols << std::endl;
	Eigen::MatrixXf mm = Eigen::MatrixXf::Zero(size, size);
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
		Eigen::MatrixXf x = ComputeMipmap(q);
#ifdef GUI_DEBUG_VERBOSE
		float xmax = x.maxCoeff();
		std::cout << "Mipmap " << mm.size() << ": min=" << x.minCoeff() << ", max=" << xmax << std::endl;
		DebugShowMatrix("x_" + boost::lexical_cast<std::string>(x.rows()), DebugDoubleMatrixSize(x,mm.size()-1), 1.0f/xmax);
#endif
		mm.push_back(x);
	}
	return mm;
}

std::vector<Eigen::Vector2f> SampleDensityImpl(const RgbdData& rgbd, const Eigen::MatrixXf& density)
{
	// compute mipmaps
	std::vector<Eigen::MatrixXf> mipmaps = ComputeMipmaps(density, 1);
	std::vector<Eigen::MatrixXf> carry_mipmaps(mipmaps.size());
	for(unsigned int i=1; i<mipmaps.size(); i++) { // HACK: skip first as we never use it
		carry_mipmaps[i] = Eigen::MatrixXf::Zero(mipmaps[i].rows(), mipmaps[i].cols());
	}
	// now create pixel seeds
	std::vector<Eigen::Vector2f> seeds;
	seeds.reserve(1000);
	SampleDensityImplRec(rgbd, seeds, mipmaps, carry_mipmaps, mipmaps.size() - 1, 0, 0);
	return seeds;
}

std::vector<Cluster> SampleClustersFromDensity(const RgbdData& rgbd, const Eigen::MatrixXf& density)
{
	std::vector<Eigen::Vector2f> points = SampleDensityImpl(rgbd, density);
	// clusters 
	std::vector<Cluster> clusters(points.size());
	for(int i=0; i<clusters.size(); i++) {
		Eigen::Vector2f px = points[i];
		int x = static_cast<int>(px.x());
		int y = static_cast<int>(px.y());
		if(!rgbd.isValid(x,y)) {
			// skip
			continue;
		}
		const Point& fp = rgbd(x,y);
		Cluster& c = clusters[i];
		c.pixel = px;
		// FIXME use small neighbourhood
		c.color = fp.color;
		c.position = fp.position;
		c.normal = fp.normal;
		c.cluster_radius_px = fp.cluster_radius_px;
		c.valid = true;
	}
	return clusters;
}

FramePtr CreateFrame(int time, const RgbdData& rgbd, const std::vector<Cluster>& clusters)
{
	constexpr float VERY_LARGE_DISTANCE = 1000000.0f;
	FramePtr p = std::make_shared<Frame>();
	p->time = time;
	p->rgbd = rgbd;
	p->clusters.resize(clusters.size());
	for(int i=0; i<clusters.size(); ++i) {
		p->clusters[i] = std::make_shared<Cluster>(clusters[i]);
		Cluster& c = *p->clusters[i];
		c.time = time;
		c.id = i;
	}
	p->assignment = assigment_type(rgbd.rows(), rgbd.cols());
	std::fill(p->assignment.begin(), p->assignment.end(), Assignment{0,VERY_LARGE_DISTANCE});
	return p;
}

inline float PointClusterDistance(int p_time, const Point& p, const Cluster& c)
{
	const float dt = static_cast<float>(std::abs(p_time-c.time)) / static_cast<float>(CLUSTER_TIME_RADIUS);
	const float dc = (p.color - c.color).squaredNorm();
	const float dx = (p.position - c.position).squaredNorm() / (CLUSTER_RADIUS*CLUSTER_RADIUS);
	return 0.5f*dt + 2.0f*dc + dx;
}

template<typename F>
void ClusterBox(const std::vector<FramePtr>& frames, F f)
{
	assert(frames.size() > 0);
	const int T = frames.size();
	const int NY = frames.front()->rgbd.cols();
	const int NX = frames.front()->rgbd.rows();
	// iterate over all frames
	for(int k=0; k<frames.size(); k++) {
		const std::vector<ClusterPtr>& frame_clusters = frames[k]->clusters;
		// iterate over clusters
		for(int i=0; i<frame_clusters.size(); i++) {
			const ClusterPtr& c = frame_clusters[i];
			// skip invalid clusters
			if(!c->valid) {
				continue;
			}
			// iterate over all pixels in box and compute distance
			for(int t=0; t<T; t++) {
				const RgbdData& rgbd = frames[t]->rgbd;
				assigment_type& assignment = frames[t]->assignment;
				int frame_time = frames[t]->time;
				// compute cluster radius
				float rpx = CLUSTER_RADIUS_MULT
					* c->cluster_radius_px
					* (1.0f + SPATIAL_TIME_INCREASE*static_cast<float>(std::abs(frame_time - c->time)) / CLUSTER_RADIUS);
				// compute cluster bounding box
				const int R = static_cast<float>(rpx + 0.5f);
				const int xc = static_cast<float>(c->pixel.x() + 0.5f);
				const int yc = static_cast<float>(c->pixel.y() + 0.5f);
				const int x1 = std::max(   0, xc - R);
				const int x2 = std::min(NX-1, xc + R);
				const int y1 = std::max(   0, yc - R);
				const int y2 = std::min(NY-1, yc + R);
				// iterate over box at time
				for(int y=y1; y<=y2; y++) {
					for(int x=x1; x<=x2; x++) {
						const Point& p = rgbd(x,y);
						// skip invalid points
						if(!p.valid) {
							continue;
						}
						// call functor
						f(c, frame_time, p, assignment(x,y));
					}
				}
			}
		}
	}
}

void UpdateClusterAssignment(const std::vector<FramePtr>& frames)
{
	ClusterBox(frames,
		[](const ClusterPtr& c, int p_time, const Point& p, Assignment& a) {
			const float d = PointClusterDistance(p_time, p, *c);
			if(d < a.distance) {
				a.distance = d;
				a.cluster = c;
			}
		});
}

struct ClusterCenterAccumulator
{
	int num;
	Eigen::Vector3f mean_color;
	Eigen::Vector3f mean_position;
	Eigen::Matrix3f mean_normal;

	ClusterCenterAccumulator()
	: num(0),
	  mean_color(Eigen::Vector3f::Zero()),
	  mean_position(Eigen::Vector3f::Zero()),
	  mean_normal(Eigen::Matrix3f::Zero())
	{}

	void add(const Point& p) {
		num ++;
		mean_color += p.color;
		mean_position += p.position;
		mean_normal += p.normal * p.normal.transpose();
	}

	Eigen::Vector3f computeNormal() const {
		// FIXME implement
		return Eigen::Vector3f(0.0f, 0.0f, -1.0f);
	}
};

void UpdateClusterCenters(Timeseries& timeseries)
{
	// init cluster center accumulators
	std::vector<std::vector<ClusterCenterAccumulator>> ccas(timeseries.frames.size());
	for(int t=0; t<ccas.size(); ++t) {
		ccas[t].resize(timeseries.frames[t]->clusters.size());
	}
	int t0 = timeseries.getBeginTime();
	// do cluster box
	ClusterBox(timeseries.frames,
		[&ccas,t0](const ClusterPtr& c, int p_time, const Point& p, Assignment& a) {
			if(a.cluster == c) {
				ccas[c->time - t0][c->id].add(p);
			}
		});
	// update
	for(int t=0; t<ccas.size(); ++t) {
		const auto& v = ccas[t];
		for(int k=0; k<v.size(); k++) {
			Cluster& c = *timeseries.frames[t]->clusters[k];
			const ClusterCenterAccumulator& cca = v[k];
			if(cca.num == 0) {
				c.valid = false;
			}
			if(!c.valid) {
				continue;
			}
			float scl = 1.0f / static_cast<float>(cca.num);
			// recompute
			c.color = scl * cca.mean_color;
			c.position = scl * cca.mean_position;
			c.normal = cca.computeNormal();
			c.pixel = CameraProject(c.position);
		}
	}
}

void UpdateClusters(int time, Timeseries& timeseries)
{
	// compute frame range for assignment update
	std::vector<FramePtr> frames = timeseries.getFrameRange(time-CLUSTER_TIME_RADIUS, time+CLUSTER_TIME_RADIUS+1);
	// iterate some times
	for (int k = 0; k < CLUSTER_ITERATIONS; ++k) {
		// update cluster assignment for frames in range
		UpdateClusterAssignment(frames);
		// update cluster centers
		UpdateClusterCenters(timeseries);
	}
}

void ContinuousSupervoxels::start(int rows, int cols)
{
	is_first_ = true;
	series_.frames.clear();
	inactive_clusters_.clear();
}

void ContinuousSupervoxels::step(const slimage::Image3ub& color, const slimage::Image1ui16& depth)
{
	constexpr float DEBUG_DENSITY_SCALE = 100.0f;

	// create rgbd data
	RgbdData rgbd = CreateRgbdData(color, depth);

	// compute frame density and sample frame clusters
	if(!is_first_) {
		DebugShowMatrix("last_density", last_density_, DEBUG_DENSITY_SCALE);
	}
	// computes frame target density
	Eigen::MatrixXf target_density = ComputeFrameDensity(rgbd);
	DebugShowMatrix("target_density", target_density, DEBUG_DENSITY_SCALE);
	// computes cluster sample density
	constexpr float LAMBDA = 1.0f - 1.0f / static_cast<float>(2*CLUSTER_TIME_RADIUS+1);
	Eigen::MatrixXf sample_density;
	if(is_first_) {
		// use target density
		sample_density = target_density;
	}
	else {
		// recent clusters provide density via last_density
		sample_density = target_density - LAMBDA*last_density_;
	}
	DebugShowMatrix("sample_density", sample_density, DEBUG_DENSITY_SCALE);
	// samples clusters from sample density
	std::vector<Cluster> new_clusters = SampleClustersFromDensity(rgbd, sample_density);
	// computes density of generated clusters
	Eigen::MatrixXf current_density = ComputeClusterDensity(rgbd.rows(), rgbd.cols(), new_clusters);
	DebugShowMatrix("current_density", current_density, DEBUG_DENSITY_SCALE);
#ifdef GUI_DEBUG_NORMAL
	slimage::gui::WaitForKeypress();
#endif
	// updates last density
	if(is_first_) {
		last_density_ = current_density;
	}
	else {
		last_density_ = LAMBDA*last_density_ + current_density;
	}

	// creates a frame and adds it to the series
	FramePtr frame = CreateFrame(series_.getEndTime(), rgbd, new_clusters);
	series_.add(frame);

	// purge old frames to limit time interval
	std::vector<ClusterPtr> purged_clusters = series_.purge(series_.getEndTime() - 2*CLUSTER_TIME_RADIUS - 1);
	inactive_clusters_.insert(inactive_clusters_.end(), purged_clusters.begin(), purged_clusters.end());

	// get current active time
	int t = std::max(series_.getBeginTime(), series_.getEndTime() - CLUSTER_TIME_RADIUS - 1);

	// Debug
	std::cout << "f=" << frame->time << ", t=" << t
			<< ", span=[" << series_.getBeginTime() << "," << series_.getEndTime() << "["
			<< ", clusters active=" << numActiveClusters() << "/inactive=" << numInactiveClusters()
			<< std::endl;

	// update clusters around active time
	UpdateClusters(t, series_);

	is_first_ = false;
}

int ContinuousSupervoxels::numActiveClusters() const
{
	return std::accumulate(series_.frames.begin(), series_.frames.end(), 0,
		[](int a, const FramePtr& f) { return a + f->clusters.size(); } );
}

int ContinuousSupervoxels::numInactiveClusters() const
{
	return inactive_clusters_.size();
}

std::vector<Cluster> ContinuousSupervoxels::getAllClusters() const
{
	std::vector<Cluster> result;
	result.reserve(numActiveClusters() + numInactiveClusters());
	for(const ClusterPtr& c : inactive_clusters_) {
		result.push_back(*c);
	}
	for(const FramePtr& f : series_.frames) {
		for(const ClusterPtr& c : f->clusters) {
			result.push_back(*c);
		}
	}
	return result;
}

void DebugWriteClusters(const std::string& fn, const std::vector<Cluster>& clusters)
{
	std::ofstream ofs(fn);
	for(const Cluster& c : clusters) {
		ofs << c.time
			<< "\t" << c.id
			<< "\t" << (c.valid ? 1 : 0)
			<< "\t" << c.cluster_radius_px
			<< "\t" << c.pixel.x() << "\t" << c.pixel.y()
			<< "\t" << c.color.x() << "\t" << c.color.y() << "\t" << c.color.z()
			<< "\t" << c.position.x() << "\t" << c.position.y() << "\t" << c.position.z()
			<< "\t" << c.normal.x() << "\t" << c.normal.y() << "\t" << c.normal.z()
			<< std::endl;
	}
}

}
