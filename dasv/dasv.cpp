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
 * - 2D arrays with user types use dasp::Array which behaves like Eigen::MatrixXf
 * - With (x,y) coordinates these correspondences should be used:
 *     width -> rows
 *     height -> cols
 *     m(x,y) (it is optimal to use x in the inner loop)
 *     a[y][x] (it is optimal to use x in the inner loop)
 */

#include "dasv.hpp"
#include <dasp/impl/Sampling.hpp>
#include <graphseg/IO.hpp>
#include <graphseg/Spectral.hpp>
#include <graphseg/Labeling.hpp>
#define DANVIL_ENABLE_BENCHMARK
#include <Danvil/Tools/Benchmark.h>
#include <Slimage/Gui.hpp>
#include <Slimage/IO.hpp>
#include <boost/format.hpp>
#include <boost/array.hpp>
#include <boost/graph/copy.hpp>
#include <random>
#include <set>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <assert.h>

//#define GUI_DEBUG_VERBOSE
#define GUI_DEBUG_NORMAL
//#define ENABLE_SAMPLING_DEBUG
//#define EVAL_COMPRESSION_ERROR

namespace dasv
{

#ifdef ENABLE_SAMPLING_DEBUG
slimage::Image3ub sampling_debug;
#endif

std::mt19937 random_engine;

constexpr float DEPTH_TO_Z = 0.001f;
constexpr float CENTER_X = 320.0f;
constexpr float CENTER_Y = 240.0f;
constexpr float PX_FOCAL = 528.0f;
constexpr float CLUSTER_RADIUS = 0.025f;
constexpr int CLUSTER_TIME_RADIUS = 0; // TR=15 -> 0.5 s
constexpr int CLUSTER_ITERATIONS = 1;
constexpr float CLUSTER_RADIUS_MULT = 1.7f;
constexpr uint16_t DEPTH_MIN = 0;
constexpr uint16_t DEPTH_MAX = 2000;
constexpr float LABELING_EDGE_MERGE_THRESHOLD = 1.0f;

constexpr float PI = 3.1415f;

boost::function<void(const std::string& tag, const slimage::Image3ub& img)> s_debug_image_display_callback
= [](const std::string& tag, const slimage::Image3ub& img) {
	slimage::gui::Show(tag, img, 1);
};

void DebugSetDisplayImageCallback(boost::function<void(const std::string& tag, const slimage::Image3ub& img)> f)
{
	s_debug_image_display_callback = f;
}

void DebugDisplayImage(const std::string& tag, const slimage::Image3ub& img)
{
	if(s_debug_image_display_callback) {
		s_debug_image_display_callback(tag, img);
	}
}

slimage::Image3ub DebugMatrixToImage(const Eigen::MatrixXf& mat, float min, float max)
{
	slimage::Image3ub img(mat.rows(), mat.cols());
	const int n = mat.size();
	for(int i=0; i<n; i++) {
		float v = mat.data()[i];
		float p = std::min(1.0f, std::max(0.0f, (v-min)/(max-min)));
		unsigned char c = static_cast<unsigned char>(255.0f*p);
		img[i] = {{c,c,c}};
	}
	return img;
}

void DebugDisplayImage(const std::string& tag, const Eigen::MatrixXf& mat, float min, float max)
{
	DebugDisplayImage(tag,
		DebugMatrixToImage(mat, min, max));
}

void DebugDisplayImage(const std::string& tag, const slimage::Image1ui16& img16, uint16_t min, uint16_t max)
{
	slimage::Image3ub img(img16.width(), img16.height());
	const int n = img16.size();
	for(int i=0; i<n; i++) {
		uint16_t v = img16[i];
		float p = std::min(1.0f, std::max(0.0f, static_cast<float>(v-min)/static_cast<float>(max-min)));
		unsigned char c = static_cast<unsigned char>(255.0f*p);
		img[i] = {{c,c,c}};
	}
	DebugDisplayImage(tag, img);
}

// Eigen::MatrixXf DebugDoubleMatrixSize(const Eigen::MatrixXf& mat, int n)
// {
// 	Eigen::MatrixXf last = mat;
// 	Eigen::MatrixXf result;
// 	for(int k=0; k<n; k++) {
// 		result = Eigen::MatrixXf(last.rows()*2, last.cols()*2);
// 		for(int i=0; i<result.cols(); i++)
// 			for(int j=0; j<result.rows(); j++)
// 				result(j,i) = last(j/2,i/2);
// 		last = result;
// 	}
// 	return last;
// }

inline float PointClusterTimeDistance(int t1, int t2)
{
	const int dt = t1 - t2;
	if(CLUSTER_TIME_RADIUS == 0) {
		return (dt == 0 ? 0.0f : 1000000.0f);
	}
	else {
		return static_cast<float>(dt*dt) / static_cast<float>(CLUSTER_TIME_RADIUS*CLUSTER_TIME_RADIUS);
	}
}

inline float ClusterClusterTimeDistance(int t1, int t2)
{
	const int dt = t1 - t2;
	if(CLUSTER_TIME_RADIUS == 0) {
		return (std::abs(dt) <= 1 ? 0.0f : 1000000.0f);
	}
	else {
		return std::max(0.0f,
			static_cast<float>(dt*dt) / static_cast<float>(4.0f*CLUSTER_TIME_RADIUS*CLUSTER_TIME_RADIUS) - 1.2f);
		// return static_cast<float>(dt*dt) / static_cast<float>(4.0f*CLUSTER_TIME_RADIUS*CLUSTER_TIME_RADIUS);
	}
}

/** Computes point to cluster distance */
inline float PointClusterDistance(int p_time, const Point& p, const Cluster& c)
{
	constexpr float wc = 1.0f;
	constexpr float ws = 1.0f;
	constexpr float wc0 = wc / (wc + ws);
	constexpr float ws0 = ws / (wc + ws);
	const float mc = 10.0f*(p.color - c.color).squaredNorm();
	const float mx = (p.position - c.position).squaredNorm() / (CLUSTER_RADIUS*CLUSTER_RADIUS);
	const float mt = PointClusterTimeDistance(p_time, c.time);
	const float ms = mx + mt;
	return wc0*mc + ws0*ms;
}

/** Computes cluster to cluster similarity */
inline float ClusterClusterSimilarity(const Cluster& a, const Cluster& b)
{
	constexpr float wc = 1.0f;
	constexpr float wx = 1.0f;
	constexpr float wt = 1.0f;
	const float mc = 3.0f*(a.color - b.color).squaredNorm();
	const float mx = std::max(0.0f,
		(a.position - b.position).squaredNorm() / (4.0f*CLUSTER_RADIUS*CLUSTER_RADIUS) - 1.2f);
	// const float mx = (a.position - b.position).squaredNorm() / (4.0f*CLUSTER_RADIUS*CLUSTER_RADIUS);
	const float mt = ClusterClusterTimeDistance(a.time, b.time);
	const float d = wc*mc + wx*mx + wt*mt;
	return std::exp(-d);
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

boost::array<unsigned char,3> ColorToImage(const Eigen::Vector3f& c) {
	return {{
		static_cast<unsigned char>(c.x() * 255.0f),
		static_cast<unsigned char>(c.y() * 255.0f),
		static_cast<unsigned char>(c.z() * 255.0f)
	}};
}

RgbdData CreateRgbdData(const Rgbd& data)
{
	const int NX = data.color.width();
	const int NY = data.color.height();
	RgbdData rgbd(NX, NY);
	for(int y=0, i=0; y<NY; y++) {
		for(int x=0; x<NX; x++, i++) {
			Point& point = rgbd(x,y);
			point.px = x;
			point.py = y;
			const uint16_t depth = data.depth[i];
			// valid
			point.is_valid = depth != 0 && DEPTH_MIN <= depth && depth <= DEPTH_MAX;
			if(point.is_valid) {
				const float z_over_f = DEPTH_TO_Z * static_cast<float>(depth) / PX_FOCAL;
				// RGB color
				const slimage::Pixel3ub& color = data.color[i];
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
			}
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
			if(p.is_valid) {
				// rho = r_px|^2 * pi / sqrt(||g||^2+1)
				// 1/sqrt(||g||^2+1) = n_z because g = -(n_x/n_z, n_y/n_z)
				// TODO n_z should always be negative so abs(n_z) should equal -n_z.
				const float A = p.cluster_radius_px * p.cluster_radius_px * PI * std::abs(p.normal.z());
				const float rho = 1.0f / A / static_cast<float>(2*CLUSTER_TIME_RADIUS+1);
				if(y==0&&x==0) std::cout << "rho " << rho << std::endl;
				density(x,y) = rho;
			}
			else {
				density(x,y) = 0.0f;
			}
		}
	}
	return density;
}

/** Iterates over a box with radius r centered around (sx,sy) */
template<typename F>
void Box(int rows, int cols, float sx, float sy, float r, F f)
{
	// box range
	const int xmin = std::max<int>(static_cast<float>(sx - r + 0.5f), 0);
	const int xmax = std::min<int>(static_cast<float>(sx + r + 0.5f), rows - 1);
	const int ymin = std::max<int>(static_cast<float>(sy - r + 0.5f), 0);
	const int ymax = std::min<int>(static_cast<float>(sy + r + 0.5f), cols - 1);
	// iterate over box pixels
	for(int yi=ymin; yi<=ymax; ++yi) {
		for(int xi=xmin; xi<=xmax; ++xi) {
			f(xi, yi);
		}
	}
}

// Eigen::MatrixXf ComputeClusterDensity(int rows, int cols, const std::vector<Cluster>& clusters)
// {
// 	// range R of kernel is s.t. phi(x) >= 0.01 * phi(0) for all x <= R
// 	const float cRange = 1.21f; // BlueNoise::KernelFunctorInverse(0.01f);
// 	Eigen::MatrixXf density = Eigen::MatrixXf::Zero(rows, cols);
// 	for(const Cluster& c : clusters) {
// 		if(!c.valid) {
// 			continue;
// 		}
// 		const float rho = 1.0f / (c.cluster_radius_px*c.cluster_radius_px*PI*std::abs(c.normal.z()));
// 		// kernel influence range
// 		const float r = cRange / std::sqrt(rho);
// 		// write kernel
// 		// seed corresponds to a kernel at position (x,y) with sigma = rho(x,y)^(-1/2)
// 		float sxf = c.pixel.x();
// 		float syf = c.pixel.y();
// 		Box(rows, cols, sxf, syf, r,
// 			[&density, sxf, syf, rho](int xi, int yi) {
// 				const float dx = static_cast<float>(xi) - sxf;
// 				const float dy = static_cast<float>(yi) - syf;
// 				const float d2 = dx*dx + dy*dy;
// 				const float delta = rho * std::exp(-PI*rho*d2);// BlueNoise::KernelFunctorSquare(rho*d2);
// 				density(xi, yi) += delta / static_cast<float>(CLUSTER_TIME_RADIUS);
// 		});
// 	}
// 	return density;
// }

Eigen::MatrixXf ComputeSeriesDensity(int ftime, const std::vector<FramePtr>& frames)
{
	constexpr float OMA = 0.618034f; // golden ratio ... ?
	constexpr float KERNEL_R_MULT = 1.21f; // ? ...
	const int rows = frames.front()->rgbd.rows();
	const int cols = frames.front()->rgbd.cols();
	Eigen::MatrixXf density = Eigen::MatrixXf::Zero(rows, cols);
	for(int i=0; i<frames.size(); i++) {
		for(const Cluster& c : frames[i]->clusters) {
			const float sx2 = PI*c.cluster_radius_px*c.cluster_radius_px*std::abs(c.normal.z());
			const float sx2_inv = 1.0f / sx2;
			const float st = static_cast<float>(2*CLUSTER_TIME_RADIUS+1);
			const float dt_over_st = static_cast<float>(ftime - frames[i]->time) / st;
			const float dt2_st2 = dt_over_st*dt_over_st;
			const float A = 1.0f / (sx2 * st);
			const float sxf = c.pixel.x();
			const float syf = c.pixel.y();
			const float kernel_r = KERNEL_R_MULT * std::sqrt(sx2);
			Box(rows, cols, sxf, syf, kernel_r,
				[&density, sxf, syf, A, sx2_inv, dt2_st2](int xi, int yi) {
					const float dx = static_cast<float>(xi) - sxf;
					const float dy = static_cast<float>(yi) - syf;
					const float d2 = dx*dx + dy*dy;
					const float delta = OMA * A * std::exp(-OMA*PI*(d2*sx2_inv + dt2_st2));
					density(xi, yi) += delta;
			});
		}
	}
	return density;
}

std::vector<Cluster> SampleClustersFromDensity(const RgbdData& rgbd, const Eigen::MatrixXf& density)
{
	std::vector<dasp::Seed> seeds = dasp::FindSeedsDepthMipmapFS(rgbd, density);
	// clusters 
	std::vector<Cluster> clusters(seeds.size());
	for(int i=0; i<clusters.size(); i++) {
		const auto& seed = seeds[i];
		const int x = seed.x;
		const int y = seed.y;
		if(!rgbd.isValid(x,y)) {
			// skip
			continue;
		}
		const Point& fp = rgbd(x,y);
		Cluster& c = clusters[i];
		c.pixel = Eigen::Vector2f(static_cast<float>(x),static_cast<float>(y));
		// FIXME use small neighbourhood
		c.color = fp.color;
		c.position = fp.position;
		c.normal = fp.normal;
		c.cluster_radius_px = fp.cluster_radius_px;
		c.valid = true;
		c.label = -1;
	}
	return clusters;
}

std::shared_ptr<ClusterContainer> ClusterList::s_storage_ = std::make_shared<ClusterContainer>();

FramePtr CreateFrame(int time, const RgbdData& rgbd, const std::vector<Cluster>& clusters)
{
	constexpr float VERY_LARGE_DISTANCE = 1000000.0f;
	FramePtr p = std::make_shared<Frame>();
	p->time = time;
	p->rgbd = rgbd;
	for(int i=0; i<clusters.size(); ++i) {
		Cluster& c = p->clusters.addCluster(clusters[i]);
		c.time = time;
	}
	p->assignment = FrameAssignment(rgbd.rows(), rgbd.cols());
	std::fill(p->assignment.begin(), p->assignment.end(), Assignment::Empty());
	return p;
}

/** Iterates over space/time box of pixels which are in range of a cluster */
template<typename F>
void ClusterBox(const std::vector<FramePtr>& frames, F f)
{
	const int T = frames.size();
	const int NY = frames.front()->rgbd.cols();
	const int NX = frames.front()->rgbd.rows();
	// iterate over all frames
	for(int k=0; k<frames.size(); k++) {
		// iterate over clusters
		for(const Cluster& c : frames[k]->clusters) {
			// skip invalid clusters
			if(!c.valid) {
				continue;
			}
			// compute cluster radius
			const float rpx = CLUSTER_RADIUS_MULT * c.cluster_radius_px;
			// iterate over all pixels in box
			for(int t=0; t<T; t++) {
				const RgbdData& rgbd = frames[t]->rgbd;
				FrameAssignment& assignment = frames[t]->assignment;
				const int frame_time = frames[t]->time;
				// iterate over cluster box
				Box(NX, NY, c.pixel.x(), c.pixel.y(), rpx,
					[&f, &rgbd, &c, &frame_time, &assignment](int x, int y) {
						const Point& p = rgbd(x,y);
						// skip invalid points
						if(!p.is_valid) {
							return;
						}
						// call functor
						f(c, frame_time, p, assignment(x,y));
					});
			}
		}
	}
}

void UpdateClusterAssignment(const std::vector<FramePtr>& frames)
{
	ClusterBox(frames,
		[](const Cluster& c, int p_time, const Point& p, Assignment& a) {
			const float d = PointClusterDistance(p_time, p, c);
			if(d < a.distance) {
				a.distance = d;
				a.cluster_id = c.cluster_id;
			}
		});
}

struct ClusterCenterAccumulator
{
	int num;
	Eigen::Vector3f mean_color;
	Eigen::Vector3f mean_position;
//	Eigen::Matrix3f mean_normal;

	ClusterCenterAccumulator()
	: num(0)
	  ,mean_color(Eigen::Vector3f::Zero())
	  ,mean_position(Eigen::Vector3f::Zero())
//	  ,mean_normal(Eigen::Matrix3f::Zero())
	{}

	void add(const Point& p) {
		num ++;
		mean_color += p.color;
		mean_position += p.position;
//		mean_normal += p.normal * p.normal.transpose();
	}

	Eigen::Vector3f computeNormal() const {
		// FIXME implement
		return Eigen::Vector3f(0.0f, 0.0f, -1.0f);
	}
};

void UpdateClusterCenters(const std::vector<FramePtr>& frames)
{
	// prepare cluster accumulators
	std::map<cluster_id_type,ClusterCenterAccumulator> ccas;
	// fill cluster accumulators
	for(int t=0; t<frames.size(); t++) {
		const auto& f = frames[t];
		const auto& fp = f->rgbd;
		const auto& fa = f->assignment;
		const int n = fa.size();
		for(int i=0; i<n; i++) {
			const auto& a = fa[i];
			// only consider pixels with a valid assignment
			if(a.hasValidCluster()) {
				// add pixel to cluster accumulator
				ccas[a.cluster_id].add(fp[i]);
			}
		}
	}
	// mark all clusters as invalid
	for(int t=0; t<frames.size(); t++) {
		for(Cluster& c : frames[t]->clusters) {
			c.valid = false;
		}
	}
	// update cluster centers
	int t0 = frames.front()->time;
	for(const auto& p : ccas) {
		cluster_id_type cid = p.first;
		Cluster& c = ClusterList::s_storage_->at(cid);
		// only update if cluster is in time range
		int k = c.time - t0;
		if(k < 0 || frames.size() <= k)
			continue;
		// compute cluster mean and update
		const ClusterCenterAccumulator& cca = p.second;
		c.valid = true;
		assert(cca.num > 0);
		float scl = 1.0f / static_cast<float>(cca.num);
		// recompute
		c.color = scl * cca.mean_color;
		c.position = scl * cca.mean_position;
		c.normal = cca.computeNormal();
		c.pixel = CameraProject(c.position);
		c.cluster_radius_px = CLUSTER_RADIUS * PX_FOCAL / c.position.z();
	}
}

void UpdateClusterEdges(const FramePtr& frame, const FramePtr& prev)
{
	const int rows = frame->rgbd.rows();
	const int cols = frame->rgbd.cols();
	const FrameAssignment& assignment = frame->assignment;
	const FrameAssignment& at1 = prev->assignment;
	std::set<Edge> edges;
//	edges.reserve(6*frame->clusters.size()); // guess number of edges
	for(int y=0; y<cols-1; y++) {
		for(int x=0; x<rows-1; x++) {
			cluster_id_type c = assignment(x,y).cluster_id;
			if(c == INVALID_CLUSTER_ID) continue;
			cluster_id_type cy1 = assignment(x,y+1).cluster_id;
			cluster_id_type cx1 = assignment(x+1,y).cluster_id;
			cluster_id_type ct1 = at1(x,y).cluster_id;
			if(cy1 != INVALID_CLUSTER_ID && c != cy1) edges.insert({c,cy1});
			if(cx1 != INVALID_CLUSTER_ID && c != cx1) edges.insert({c,cx1});
			if(ct1 != INVALID_CLUSTER_ID && c != ct1) edges.insert({c,ct1});
		}
	}
	// return as vector
	frame->edges = std::vector<Edge>(edges.begin(), edges.end());
}

void UpdateClusterEdges(const FramePtr& frame)
{
	UpdateClusterEdges(frame, frame);
}

void UpdateClusters(const std::vector<FramePtr>& frames)
{
	// iterate some times
	for (int k = 0; k < CLUSTER_ITERATIONS; ++k) {
		// update cluster assignment for frames in range
		UpdateClusterAssignment(frames);
		// update cluster centers
		UpdateClusterCenters(frames);
	}
}

ClusterGraph ComputeClusterGraphImpl(const std::set<int>& cluster_ids, const std::set<Edge>& edges)
{
	// create graph
	ClusterGraph G(cluster_ids.size());

	// create vertices
	int i = 0;
	std::map<int,ClusterGraph::vertex_descriptor> cid_to_vid;
	for(int cid : cluster_ids) {
		G[i] = ClusterList::s_storage_->at(cid);
		cid_to_vid[cid] = i;
		i++;
	}

	// create edges
	for(const Edge& e : edges) {
		auto ea_it = cid_to_vid.find(e.a);
		auto eb_it = cid_to_vid.find(e.b);
		ClusterGraph::edge_descriptor eid;
		bool ok;
		boost::tie(eid,ok) = boost::add_edge(ea_it->second, eb_it->second, G);
		assert(ok);
		boost::put(boost::edge_weight, G, eid, 1.0f);
	}

	return G;
}


ClusterGraph ComputeClusterGraphStrict(const std::vector<FramePtr>& frames)
{
	// compute edge set
	std::set<Edge> edges;
	for(const FramePtr& f : frames) {
		edges.insert(f->edges.begin(), f->edges.end());
	}

	// collect cluster ids
	std::set<int> cluster_ids;
	for(const FramePtr& f : frames) {
		for(const Cluster& c : f->clusters) {
			cluster_ids.insert(c.cluster_id);
		}
	}

	return ComputeClusterGraphImpl(cluster_ids, edges);
}

ClusterGraph ComputeClusterGraph(const std::vector<FramePtr>& frames)
{
	// compute edge set
	std::set<Edge> edges;
	for(const FramePtr& f : frames) {
		edges.insert(f->edges.begin(), f->edges.end());
	}

	// collect cluster ids
	std::set<int> cluster_ids;
	for(const Edge& e : edges) {
		cluster_ids.insert(e.a);
		cluster_ids.insert(e.b);
	}

	return ComputeClusterGraphImpl(cluster_ids, edges);
}

void ComputeClusterGraphEdgeWeights(ClusterGraph& graph)
{
	for(auto eid : as_range(boost::edges(graph))) {
		auto ea = boost::source(eid, graph);
		auto eb = boost::target(eid, graph);
		const Cluster& ca = graph[ea];
		const Cluster& cb = graph[eb];
		const float sim = ClusterClusterSimilarity(ca, cb);
		boost::put(boost::edge_weight, graph, eid, sim);
	}
}

void IOWriteClusters(const std::string& fn, const std::vector<Cluster>& clusters)
{
	std::ofstream ofs(fn);
	for(const Cluster& cp : clusters) {
		ofs << cp << std::endl;
	}
}

std::vector<Cluster> IOReadClusters(const std::string& fn)
{
	std::vector<Cluster> clusters;
	std::ifstream ifs(fn);
	while(ifs) {
		Cluster c;
		ifs
			>> c.time >> c.cluster_id >> c.valid >> c.cluster_radius_px
			>> c.pixel.x() >> c.pixel.y()
			>> c.color.x() >> c.color.y() >> c.color.z()
			>> c.position.x() >> c.position.y() >> c.position.z()
			>> c.normal.x() >> c.normal.y() >> c.normal.z();
		clusters.push_back(c);
	}
	return clusters;
}

void IOWriteEdges(const std::string& fn, const std::vector<Edge>& edges)
{
	std::ofstream ofs(fn);
	for(const Edge& e : edges) {
		ofs << e.a << "\t" << e.b << std::endl;
	}
}

void IOWriteGraph(const std::string& fn_vertices, const std::string& fn_edges, const ClusterGraph& graph)
{
	// write clusters
	std::ofstream ofsv(fn_vertices);
	for(auto vid : as_range(boost::vertices(graph))) {
		ofsv << graph[vid] << std::endl;
	}
	// write edges
	graphseg::WriteEdges(fn_edges, graph, boost::get(boost::edge_weight, graph));
}

ClusterGraph IOReadGraph(const std::string& fn_vertices, const std::string& fn_edges)
{
	std::vector<Cluster> clusters = IOReadClusters(fn_vertices);
	ClusterGraph graph(clusters.size());
	for(std::size_t i=0; i<clusters.size(); i++) {
		graph[i] = clusters[i];
	}
	graphseg::ReadEdges(fn_edges, graph, boost::get(boost::edge_weight, graph));
	return graph;
}

void ContinuousSupervoxels::start()
{
	is_first_ = true;
	series_.frames.clear();
}

void ContinuousSupervoxels::step(const Rgbd& data)
{
	constexpr float DEBUG_DENSITY_MAX = 1.0f/800.0f;

	DANVIL_BENCHMARK_START(dasv_debug)
	DebugDisplayImage("color", data.color);
	DebugDisplayImage("depth", data.depth, 500, 3000);
	DANVIL_BENCHMARK_STOP(dasv_debug)

	// *** PREPARATIONS

	int time_new = series_.frames.empty() ? 0 : (series_.frames.back()->time + 1);
	// density range (new frame not yes added to series)
	int time_density_a = std::max(time_new - 2*CLUSTER_TIME_RADIUS, 0);
	int time_density_b = time_new - 1;
	// clustering range (new frame added to series)
	int time_clustering_a = std::max(time_new - 2*CLUSTER_TIME_RADIUS, 0);
	int time_clustering_b = time_new;
	int time_edge_update = time_new - 2*CLUSTER_TIME_RADIUS;
	bool do_edge_update = (time_edge_update >= 0);
	// labeling range
	int time_labeling_a = std::max(time_new - 4*CLUSTER_TIME_RADIUS - 1, 0);
	int time_labeling_update = time_new - 3*CLUSTER_TIME_RADIUS;
	int time_labeling_b = time_new - 2*CLUSTER_TIME_RADIUS;
	bool do_labeling = (time_labeling_update >= 0);
	// visual
	int time_visual = time_new - 4*CLUSTER_TIME_RADIUS - 1;
	bool do_visual = (time_visual >= 0);
	// purge range
	int time_purged = time_new - 4*CLUSTER_TIME_RADIUS - 1;
	bool do_purge = (time_purged >= 0);

	// Debug
	std::cout << "new=" << time_new
		<< ", density=[" << time_density_a << "," << time_density_b << "]"
		<< ", clustering=[" << time_clustering_a << "," << time_clustering_b << "]"
		<< ", edge=" << time_edge_update
		<< ", labeling=[" << time_labeling_a << "," << time_labeling_b << "]"
		<< ", labeling up=" << time_labeling_update
		<< ", visualization=" << time_visual
		<< ", purged=" << time_purged
//		<< ", clusters active=" << numActiveClusters() << "/inactive=" << numInactiveClusters()
		<< std::endl;

	// *** ADD NEW FRAME

	// create rgbd data
	DANVIL_BENCHMARK_START(dasv_rgbd)
	RgbdData rgbd = CreateRgbdData(data);
	DANVIL_BENCHMARK_STOP(dasv_rgbd)

	// computes frame target density
	DANVIL_BENCHMARK_START(dasv_density)
	Eigen::MatrixXf target_density = ComputeFrameDensity(rgbd);
	DANVIL_BENCHMARK_STOP(dasv_density)

	// density from all clusters up to now
	DANVIL_BENCHMARK_START(dasv_series_density)
	if(is_first_ || time_density_a > time_density_b) {
		last_density_ = Eigen::MatrixXf::Zero(rgbd.rows(), rgbd.cols());
	}
	else {
		std::vector<FramePtr> density_frames = series_.getFrameRange(time_density_a, time_density_b);
		last_density_ = ComputeSeriesDensity(time_new, density_frames);
	}
	DANVIL_BENCHMARK_STOP(dasv_series_density)

	// compute sample density
	DANVIL_BENCHMARK_START(dasv_sampling)
	// -> avoid sampling at same positions as last frame!
	// -> but total density shall not change!
	float target_density_sum = target_density.sum();
	float last_density_sum = last_density_.sum();
	float mult = 1.0f + last_density_sum / target_density_sum;
	Eigen::MatrixXf sample_density = mult*target_density - last_density_;
	// samples clusters from sample density
	std::vector<Cluster> new_clusters = SampleClustersFromDensity(rgbd, sample_density);
	DANVIL_BENCHMARK_STOP(dasv_sampling)

#ifdef ENABLE_SAMPLING_DEBUG
	DebugDisplayImage("sampling debug", sampling_debug);
	std::cout << "sample_density_sum=" << sample_density.sum() << std::endl;
#endif

	std::cout << "Num clusters: " << new_clusters.size() << std::endl;
	// // computes density of generated clusters
	// Eigen::MatrixXf current_density = ComputeClusterDensity(rgbd.rows(), rgbd.cols(), new_clusters);
	// debug
#ifdef ENABLE_SAMPLING_DEBUG
	//if(!is_first_) {
		DebugDisplayImage("last_density", last_density_, 0.0f, DEBUG_DENSITY_MAX);
		std::cout << "last_density_sum=" << last_density_.sum() << std::endl;
	//}
	DebugDisplayImage("target_density", target_density, 0.0f, DEBUG_DENSITY_MAX);
	std::cout << "target_density_sum=" << target_density.sum() << std::endl;
	// DebugShowMatrix("current_density", current_density, 0.0f, DEBUG_DENSITY_MAX);
	// std::cout << "current_density_sum=" << current_density.sum() << std::endl;
#endif

	// creates a frame and adds it to the series
	DANVIL_BENCHMARK_START(dasv_create_frame)
	FramePtr new_frame = CreateFrame(time_new, rgbd, new_clusters);
	series_.add(new_frame);
	DANVIL_BENCHMARK_STOP(dasv_create_frame)

	// *** CLUSTERING

#ifdef ENABLE_SAMPLING_DEBUG
	{
		// compute frame density of clustering range
		DANVIL_BENCHMARK_START(dasv_series_density)
		Eigen::MatrixXf density_now = ComputeSeriesDensity(
			timeseries.getFrameRange(time_clustering_a, time_clustering_b);
		DANVIL_BENCHMARK_STOP(dasv_series_density)
		DebugDisplayImage("density_now", density_now, 0.0f, DEBUG_DENSITY_MAX);
		std::cout << "density_now_sum=" << density_now.sum() << std::endl;
	}
#endif

	// update pixel to cluster assignement
	DANVIL_BENCHMARK_START(dasv_update_clusters)
	std::vector<FramePtr> clustering_frames = series_.getFrameRange(time_clustering_a, time_clustering_b);
	UpdateClusters(clustering_frames);
	DANVIL_BENCHMARK_STOP(dasv_update_clusters)

	// compute edges for last frame of clustering range
	if(do_edge_update) {
		FramePtr f_edge_up = series_.getFrame(time_edge_update);
		if(time_edge_update == 0) {
			// compute inner edges (no previous frame)
			UpdateClusterEdges(f_edge_up);
		}
		else {
			// compute inner edges and edges to previous frame
			UpdateClusterEdges(f_edge_up, series_.getFrame(time_edge_update-1));
		}
		std::cout << "num new edges=" << f_edge_up->edges.size() << std::endl;
	}

	// *** LABELING

	if(do_labeling) {
	 	DANVIL_BENCHMARK_START(dasv_labeling)

		// compute graph for labeling frame range
		std::vector<FramePtr> labeling_frames = series_.getFrameRange(time_labeling_a, time_labeling_b);
		ClusterGraph graph = ComputeClusterGraph(labeling_frames);

		// compute cluster similarity
		ComputeClusterGraphEdgeWeights(graph);

		// store graph
		graph_ = graph;
		std::cout << "Cluster Graph: num_vertices=" << boost::num_vertices(graph_)
			<< ", num_edges=" << boost::num_edges(graph_)
			<< std::endl;

		// apply graph segmentation
		graphseg::SpectralGraph spectral;
		boost::copy_graph(graph, spectral,
				boost::edge_copy(
					[&graph,&spectral](ClusterGraph::edge_descriptor src, typename graphseg::SpectralGraph::edge_descriptor dst) {
						boost::put(boost::edge_weight, spectral, dst,
							boost::get(boost::edge_weight, graph, src));
					}
		));
		graphseg::SpectralGraph solved = graphseg::SolveSpectral(spectral, 12);

		graph_segmented_ = ClusterGraph();
		boost::copy_graph(solved, graph_segmented_,
				boost::vertex_copy(
					[&graph_,&graph_segmented_](typename graphseg::SpectralGraph::vertex_descriptor src, ClusterGraph::vertex_descriptor dst) {
						graph_segmented_[dst] = graph_[dst];
					}
				)
				.edge_copy(
					[&solved,&graph_segmented_](typename graphseg::SpectralGraph::edge_descriptor src, ClusterGraph::edge_descriptor dst) {
						boost::put(boost::edge_weight, graph_segmented_, dst,
							boost::get(boost::edge_weight, solved, src));
					}
		));

		// compute labels (supervised)
		std::vector<int> labeling = graphseg::ComputeSegmentLabels_UCM_Supervised(
			solved,
			boost::get(boost::edge_weight, solved),
			boost::get(&Cluster::label, graph),
			LABELING_EDGE_MERGE_THRESHOLD);

		// write back labels to clusters
		FramePtr label_update_frame = series_.getFrame(time_labeling_update);
		for(auto vid : as_range(boost::vertices(graph))) {
			Cluster& c = ClusterList::s_storage_->at(graph[vid].cluster_id);
			if(c.time == time_labeling_update) {
				c.label = labeling[vid];
			}
		}

	 	DANVIL_BENCHMARK_STOP(dasv_labeling)
	}

// 	DANVIL_BENCHMARK_START(dasv_graph)
// 	// purge old frames to limit time interval
// 	std::vector<FramePtr> purged_frames = series_.purge(series_.getEndTime() - 2*CLUSTER_TIME_RADIUS - 1);
// 	// store purged frames in graph
// 	{
// 		// count number of clusters
// 		int num_clusters = ClusterList::s_storage_->data().size();
// 		// create vertices
// 		ClusterGraph G(num_clusters);
// //		std::map<int,ClusterGraph::vertex_descriptor> cid_to_vid;
// 		int i = 0;
// 		for(const Cluster& c : ClusterList::s_storage_->data()) {
// 			G[i] = c;
// 			if(c.cluster_id != i) {
// 				std::cerr << "ERROR: Cluster ID does not match array position!" << std::endl;
// 			}
// //			cid_to_vid[c.cluster_id] = i;
// 			i++;
// 		}
// 		// create edges
// 		for(auto eid : as_range(boost::edges(graph_))) {
// 			ClusterGraph::edge_descriptor neid;
// 			bool ok;
// 			boost::tie(neid,ok) = boost::add_edge(
// 				boost::source(eid,graph_),
// 				boost::target(eid,graph_),
// 				G);
// 			boost::put(boost::edge_weight, G, neid,
// 				boost::get(boost::edge_weight, graph_, eid));
// 		}
// 		std::vector<Edge> edges;
// 		for(const FramePtr& f : purged_frames) {
// 			edges.insert(edges.end(), f->edges.begin(), f->edges.end());
// 		}
// 		for(const Edge& e : edges) {
// //			auto ea_it = cid_to_vid.find(e.a);
// //			auto eb_it = cid_to_vid.find(e.b);
// //			if(ea_it == cid_to_vid.end() || eb_it == cid_to_vid.end()) {
// //				std::cerr << "ERROR: Invalid edge!" << std::endl;
// //			}
// 			ClusterGraph::edge_descriptor eid;
// 			bool ok;
// //			boost::tie(eid,ok) = boost::add_edge(ea_it->second, eb_it->second, G);
// 			boost::tie(eid,ok) = boost::add_edge(e.a, e.b, G);
// 			if(ok) {
// 				boost::put(boost::edge_weight, G, eid, 1.0f);
// 			}
// 		}
// 		ComputeClusterGraphEdgeWeights(G);
// 		// ready
// 		graph_ = G;
// 		std::cout << "Graph: num_vertices=" << boost::num_vertices(graph_)
// 				<< ", num_edges=" << boost::num_edges(graph_)
// 				<< std::endl;
// 	}
// 	DANVIL_BENCHMARK_STOP(dasv_graph)

	// *** VISUALIZATION

#ifdef GUI_DEBUG_NORMAL
	{
		// superpixel image
		DANVIL_BENCHMARK_START(dasv_debug)

		if(do_visual) {
			FramePtr frame = series_.getFrame(time_visual);

#ifdef EVAL_COMPRESSION_ERROR
		{
			// computes compression error
			Eigen::Vector2f compression_error = EvaluateComputeCompressionError(frame);
			Eigen::Vector2f ref_compression_error = EvaluateComputeDownsampleCompressionError(frame);
			std::cout << "Compression Error (t=" << time_clustering_vis << "): " << compression_error.transpose() << "(ref=" << ref_compression_error.transpose() << ")" << std::endl;
		}
#endif

			boost::format fmt_col("/tmp/dasv/%05d_color.png");
			slimage::Image3ub img_col = DebugPlotClusters(frame, {PlotStyle::Color, PlotStyle::ClusterBorder});
			DebugDisplayImage("sv color", img_col);
			slimage::Save(img_col, (fmt_col % frame->time).str());

			// boost::format fmt_age("/tmp/dasv/%05d_age.png");
			// slimage::Image3ub img_age = DebugPlotClusters(frame, {PlotStyle::Age, PlotStyle::ClusterBorder});
			// DebugDisplayImage("sv age", img_age);
			// slimage::Save(img_age, (fmt_age % frame->time).str());

			// slimage::Image3ub img_adist = DebugPlotClusters(frame, {PlotStyle::AssignmentDistance});
			// DebugDisplayImage("sv assignment dist", img_adist);

			boost::format fmt_label("/tmp/dasv/%05d_label.png");
			slimage::Image3ub img_label = DebugPlotClusters(frame, {PlotStyle::Label});
			DebugDisplayImage("sv label", img_label);
			slimage::Save(img_label, (fmt_label % frame->time).str());
		}

		// cluster graph
		boost::format fmt_clusters("/tmp/dasv/%05d_clusters.tsv");
		boost::format fmt_edges("/tmp/dasv/%05d_edges.tsv");
		IOWriteGraph(
			(fmt_clusters % time_new).str(),
			(fmt_edges % time_new).str(),
			graph_);

		// cluster graph
		boost::format fmt_clusters2("/tmp/dasv/%05d_clusters2.tsv");
		boost::format fmt_edges2("/tmp/dasv/%05d_edges2.tsv");
		IOWriteGraph(
			(fmt_clusters2 % time_new).str(),
			(fmt_edges2 % time_new).str(),
			graph_segmented_);

		DANVIL_BENCHMARK_STOP(dasv_debug)
	}
#endif

	// *** PURGE UNNEEDED FRAMES

	if(do_purge) {
		// remove first frame
		series_.frames.erase(series_.frames.begin(), series_.frames.begin()+1);
	}

	// *** FINISHED

#ifdef GUI_DEBUG_NORMAL
	if((time_new+1) % 10 == 0) {
		DANVIL_BENCHMARK_PRINTALL_COUT
	}
#endif

	is_first_ = false;
}

int ContinuousSupervoxels::numActiveClusters() const
{
	return std::accumulate(series_.frames.begin(), series_.frames.end(), 0,
		[](int a, const FramePtr& f) { return a + f->clusters.size(); } );
}

int ContinuousSupervoxels::numInactiveClusters() const
{
	return ClusterList::s_storage_->data().size() - numActiveClusters();
}

int ContinuousSupervoxels::numClusters() const
{
	return numActiveClusters() + numInactiveClusters();
}

std::vector<Cluster> ContinuousSupervoxels::getAllClusters() const
{
	return ClusterList::s_storage_->data();
}

template<typename ColorFunc>
void DebugPlotSuperpixels(slimage::Image3ub& img, const FramePtr& frame, ColorFunc cf)
{
	const FrameAssignment& assignment = frame->assignment;
	const int rows = frame->rgbd.rows();
	const int cols = frame->rgbd.cols();
	img.resize(rows, cols);
	for(int y=0; y<cols; y++) {
		for(int x=0; x<rows; x++) {
			slimage::Pixel3ub color;
			const auto& a = assignment(x,y);
			if(a.hasValidCluster()) {
				color = cf(a);
			}
			else {
				color = (x%2==y%2)
					? slimage::Pixel3ub{{96,0,96}}
					: slimage::Pixel3ub{{0,0,0}};
			}
			img(x,y) = color;
		}
	}
}

void DebugPlotColor(slimage::Image3ub& img, const FramePtr& frame)
{
	DebugPlotSuperpixels(img, frame,
		[](const Assignment& a) -> slimage::Pixel3ub {
			const Cluster& c = a.getCluster();
			// cluster color for pixel
			const auto pc = ColorToImage(c.color);
			return {{pc[0],pc[1],pc[2]}};
		});
}

void DebugPlotAge(slimage::Image3ub& img, const FramePtr& frame)
{
	const auto ftime = frame->time;
	DebugPlotSuperpixels(img, frame,
		[ftime](const Assignment& a) -> slimage::Pixel3ub {
			const Cluster& c = a.getCluster();
			const int dt = ftime - c.time;
			const int q = 
				(CLUSTER_TIME_RADIUS == 0)
				? (
					(dt == 0 ? 0 : (dt < 0 ? -512 : +512))
				)
				: (
					(dt*255)/CLUSTER_TIME_RADIUS
				);
			if(q < -255) {
				return {{ 0,96,0 }};
			}
			else if(q > +255) {
				//return {{ (unsigned char)((510-q)/2), 0, 0 }};
				return {{ 96,0,0 }};
			}
			else if(q < 0) {
				return {{ (unsigned char)(255+q), 255, 0 }};
			}
			else {
				return {{ 255, (unsigned char)(255-q), 0 }};
			}
		});
}

void DebugPlotAssignmentDistance(slimage::Image3ub& img, const FramePtr& frame)
{
	DebugPlotSuperpixels(img, frame,
		[](const Assignment& a) -> slimage::Pixel3ub {
			const float x = std::min(+1.0f, std::max(-1.0f, a.distance/0.5f - 1.0f));
			const int q = static_cast<int>(x*255.0f);
			if(q < 0) {
				return {{ 0, 0, (unsigned char)(-q) }};
			}
			else {
				return {{ (unsigned char)(q), 0, 0 }};
			}
		});
}

void DebugPlotLabel(slimage::Image3ub& img, const FramePtr& frame)
{
	static std::map<int, slimage::Pixel3ub> label_colors;
	label_colors[-1] = {{255,0,255}};
	label_colors[0] = {{64,64,192}};
	DebugPlotSuperpixels(img, frame,
		[&label_colors](const Assignment& a) -> slimage::Pixel3ub {
			const int label = a.getCluster().label;
			// find color
			slimage::Pixel3ub color;
			auto it = label_colors.find(label);
			if(it == label_colors.end()) {
				std::uniform_int_distribution<int> dist(0,255);
				color = {{
					static_cast<unsigned char>(dist(random_engine)),
					static_cast<unsigned char>(dist(random_engine)),
					static_cast<unsigned char>(dist(random_engine)) }};
				// create new random color
				label_colors[label] = color;
			}
			else {
				color = it->second;
			}
			// cluster color for pixel
			return color;
		});
}

void DebugPlotClusterBorder(slimage::Image3ub& img, const FramePtr& frame)
{
	const FrameAssignment& assignment = frame->assignment;
	const int rows = frame->rgbd.rows();
	const int cols = frame->rgbd.cols();
	img.resize(rows, cols);
	for(int y=1; y<cols-1; y++) {
		for(int x=1; x<rows-1; x++) {
			cluster_id_type cid = assignment(x,y).cluster_id;
			if(    cid != assignment(x,y-1).cluster_id
				|| cid != assignment(x-1,y).cluster_id
				|| cid != assignment(x,y+1).cluster_id
				|| cid != assignment(x+1,y).cluster_id
			) {
				const slimage::Pixel3ub& v = img(x,y);
				unsigned char cr = 255 - v[0];
				unsigned char cg = 255 - v[1];
				unsigned char cb = 255 - v[2];
				img(x,y) = {{cr, cg, cb}};
			}
		}
	}
}

void DebugPlotClusters(slimage::Image3ub& img, const FramePtr& frame, PlotStyle style)
{
	switch(style) {
		case PlotStyle::Color:
			return DebugPlotColor(img, frame);
		case PlotStyle::Age:
			return DebugPlotAge(img, frame);
		case PlotStyle::AssignmentDistance:
			return DebugPlotAssignmentDistance(img, frame);
		case PlotStyle::Label:
			return DebugPlotLabel(img, frame);
		case PlotStyle::ClusterBorder:
			return DebugPlotClusterBorder(img, frame);
	}
}

slimage::Image3ub DebugPlotClusters(const FramePtr& frame, const std::vector<PlotStyle>& styles)
{
	slimage::Image3ub img;
	for(PlotStyle ps : styles) {
		DebugPlotClusters(img, frame, ps);
	}
	return img;
}


Eigen::Vector2f EvaluateComputeCompressionError(const FramePtr& frame)
{
	const int n = frame->rgbd.size();
	const RgbdData& rgbd = frame->rgbd;
	const FrameAssignment& assignment = frame->assignment;
	// compute mean of all pixels
	Eigen::Vector3f pixel_mean_color = Eigen::Vector3f::Zero();
	Eigen::Vector3f pixel_mean_position = Eigen::Vector3f::Zero();
	int num_pixels = 0;
	for(int i=0; i<n; i++) {
		const Point& p = rgbd[i];
		if(!p.is_valid) continue;
		pixel_mean_color += p.color;
		pixel_mean_position += p.position;
		num_pixels ++;
	}
	assert(num_pixels > 0);
	pixel_mean_color /= static_cast<float>(num_pixels);
	pixel_mean_position /= static_cast<float>(num_pixels);
	// compute errors
	float cluster_error_color = 0.0f;
	float cluster_error_position = 0.0f;
	float pixel_error_color = 0.0f;
	float pixel_error_position = 0.0f;
	for(int i=0; i<n; i++) {
		const Assignment& a = assignment[i];
		if(!a.hasValidCluster())
			continue;
		const Point& p = rgbd[i];
		if(!p.is_valid)
			continue;
		const Cluster& c = a.getCluster();
		cluster_error_color += (c.color - pixel_mean_color).squaredNorm();
		cluster_error_position += (c.position - pixel_mean_position).squaredNorm();
		pixel_error_color += (p.color - pixel_mean_color).squaredNorm();
		pixel_error_position += (p.position - pixel_mean_position).squaredNorm();
	}
	return {
		cluster_error_color / pixel_error_color,
		cluster_error_position / pixel_error_position
	};
}

Eigen::Vector2f EvaluateComputeDownsampleCompressionError(const FramePtr& frame)
{
	const RgbdData& rgbd = frame->rgbd;
	const FrameAssignment& assignment = frame->assignment;
	// compute mean of all pixels
	Eigen::Vector3f pixel_mean_color = Eigen::Vector3f::Zero();
	Eigen::Vector3f pixel_mean_position = Eigen::Vector3f::Zero();
	int num_pixels = 0;
	const int n = rgbd.size();
	for(int i=0; i<n; i++) {
		const Point& p = rgbd[i];
		if(!p.is_valid) continue;
		pixel_mean_color += p.color;
		pixel_mean_position += p.position;
		num_pixels ++;
	}
	// downsample assignment
	assert(num_pixels > 0);
	pixel_mean_color /= static_cast<float>(num_pixels);
	pixel_mean_position /= static_cast<float>(num_pixels);
	const int num_clusters = frame->clusters.size();
	const int rows = rgbd.rows();
	const int cols = rgbd.cols();
	const int sclrows = static_cast<int>(3.464f*std::sqrt(num_clusters) + 0.5f);
	const int sclcols = static_cast<int>(2.598f*std::sqrt(num_clusters) + 0.5f);
	std::cout << sclrows << " " << sclcols << std::endl;
	dasp::Array<Eigen::Vector3f> cluster_color(sclrows, sclcols, Eigen::Vector3f::Zero());
	dasp::Array<Eigen::Vector3f> cluster_position(sclrows, sclcols, Eigen::Vector3f::Zero());
	dasp::Array<int> num(sclrows, sclcols, 0);
	for(int i=0; i<cols; i++) {
		for(int j=0; j<rows; j++) {
			const Point& p = rgbd(j,i);
			if(!p.is_valid) continue;
			const int si = (i * sclcols) / cols;
			const int sj = (j * sclrows) / rows;
			cluster_color(sj,si) += p.color;
			cluster_position(sj,si) += p.position;
			num(sj,si) ++;
		}
	}
	const int scln = cluster_color.size();
	for(int i=0; i<scln; i++) {
		if(num[i] == 0) continue;
		cluster_color[i] /= static_cast<float>(num[i]);
		cluster_position[i] /= static_cast<float>(num[i]);
	}
	// compute errors
	float cluster_error_color = 0.0f;
	float cluster_error_position = 0.0f;
	float pixel_error_color = 0.0f;
	float pixel_error_position = 0.0f;
	for(int i=0; i<cols; i++) {
		for(int j=0; j<rows; j++) {
			const int si = (i * sclcols) / cols;
			const int sj = (j * sclrows) / rows;
			const Point& p = rgbd(j,i);
			if(num(sj,si) == 0 || !p.is_valid) continue;
			cluster_error_color += (cluster_color(sj,si) - pixel_mean_color).squaredNorm();
			cluster_error_position += (cluster_position(sj,si) - pixel_mean_position).squaredNorm();
			pixel_error_color += (p.color - pixel_mean_color).squaredNorm();
			pixel_error_position += (p.position - pixel_mean_position).squaredNorm();
		}
	}
	return {
		cluster_error_color / pixel_error_color,
		cluster_error_position / pixel_error_position
	};
}

}
