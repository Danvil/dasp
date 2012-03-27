/*
 * Superpixel.cpp
 *
 *  Created on: Jan 26, 2012
 *      Author: david
 */

#include "Superpixels.hpp"
#include "tools/Graph.hpp"
#include "tools/RepairDepth.hpp"
#define DANVIL_ENABLE_BENCHMARK
#include <Danvil/Tools/Benchmark.h>
#include <Danvil/Tools/MoreMath.h>
#include <Danvil/Color.h>
#include <Danvil/Color/LAB.h>
#include <Danvil/Color/HSV.h>
#include <eigen3/Eigen/Eigenvalues>
#include <boost/math/constants/constants.hpp>
#include <queue>
#include <set>

//#define CLUSTER_UPDATE_MULTITHREADING
#define CREATE_DEBUG_IMAGES

//------------------------------------------------------------------------------
namespace dasp {
//------------------------------------------------------------------------------

std::map<std::string,slimage::ImagePtr> sDebugImages;

Parameters::Parameters()
{
	color_space = ColorSpaces::RGB;
	weight_color = 2.0f;
	weight_spatial = 1.0f;
	weight_normal = 3.0f;
	weight_depth = 0.0f;
	weight_image = 0.0f;
	iterations = 5;
	coverage = 1.7f;
	base_radius = 0.02f;
	count = 0;
	seed_mode = SeedModes::DepthMipmap;
	gradient_adaptive_density = true;
	is_conquer_enclaves = true;
	segment_threshold = 1.0f;
	is_repair_depth = true;
	is_smooth_depth = false;
	is_improve_seeds = false;
}

void Cluster::UpdateCenter(const ImagePoints& points, const Parameters& opt)
{
	assert(hasPoints());

	// compute mean
	Eigen::Vector3f mean_color = Eigen::Vector3f::Zero();
	Eigen::Vector3f mean_world = Eigen::Vector3f::Zero();

	for(unsigned int i : pixel_ids) {
		const Point& p = points[i];
		assert(p.isValid());
		mean_color += p.color;
		mean_world += p.world;
	}

	center.color = mean_color / float(pixel_ids.size());
	center.world = mean_world / float(pixel_ids.size());
	center.pos = opt.camera.project(center.world);
	center.depth_i16 = opt.camera.depth(center.world);

	// FIXME change or not change? (SLIC mode does not allow change!)
	//center.image_super_radius = opt.base_radius * opt.camera.scala(center.depth_i16);

	cov = PointCovariance(pixel_ids, [this,&points](unsigned int i) { return points[i].world - center.world; });
	Eigen::SelfAdjointEigenSolver<Eigen::Matrix3f> solver;
	solver.computeDirect(cov);
	ew = solver.eigenvalues();
	ev = solver.eigenvectors();

	Eigen::Vector3f normal = ev.col(0);
	if(normal[2] == 0.0f) {
		normal[2] = -0.01f;
	}
	else if(normal[2] > 0.0f) {
		normal *= -1.0f;
	}
	normal.normalize();
	center.gradient[0] = - normal[0] / normal[2];
	center.gradient[1] = - normal[1] / normal[2];
	center.circularity = std::abs(normal[2]); // == 1.0f / std::sqrt(center.gradient.squaredNorm() + 1.0f);

//	center.normal = points(center.pos).normal;
//	center.gradient = points(center.pos).gradient;

	// eigenvalues are the square of the standard deviation!

	thickness = cSigmaScale * std::sqrt(std::abs(ew(0))) * 2.0f;

	circularity = std::sqrt(std::abs(ew(1) / ew(2)));

	eccentricity = std::sqrt(1.0f - std::abs(ew(1) / ew(2)));

	area_quotient = cSigmaScale * cSigmaScale * std::sqrt(std::abs(ew(1) * ew(2))) / (opt.base_radius * opt.base_radius);

}

void Cluster::ComputeExt(const ImagePoints& points, const Parameters& opt)
{
	std::vector<float> dist;
	dist.reserve(pixel_ids.size());
	for(unsigned int id : pixel_ids) {
		Eigen::Vector3f p = points[id].world - center.world;
		float d = p.dot(center.computeNormal());
		dist.push_back(d);
	}
	unsigned int nth = static_cast<unsigned int>(cPercentage * static_cast<float>(dist.size()));
	std::nth_element(dist.begin(), dist.begin() + nth, dist.end());
	thickness_plane = (*(dist.begin() + nth));

	{
		area_actual = 0.0f;
		area_expected = 0.0f;
		float error_area = 0;
		int cx = center.spatial_x();
		int cy = center.spatial_y();
		int R = int(center.image_super_radius * opt.coverage * 1.5f);
		unsigned int xmin = std::max(0, cx - R);
		unsigned int xmax = std::min(int(points.width()-1), cx + R);
		unsigned int ymin = std::max(0, cy - R);
		unsigned int ymax = std::min(int(points.height()-1), cy + R);
		for(unsigned int y=ymin; y<=ymax; y++) {
			for(unsigned int x=xmin; x<=xmax; x++) {
				unsigned int i = points.index(x,y);
				const Point& p = points[i];
				Eigen::Vector3f t = p.world - center.world;
				float lam = t.dot(center.computeNormal());
				float d2 = t.squaredNorm() - lam*lam;
				bool expected = lam < thickness_plane && d2 < opt.base_radius*opt.base_radius;
				bool actual = std::find(pixel_ids.begin(), pixel_ids.end(), i) != pixel_ids.end();
				float size_of_a_px = p.depth() / opt.camera.focal;
				float area_of_a_px = size_of_a_px*size_of_a_px / p.circularity;
				if(actual) {
					area_actual += area_of_a_px;
				}
				if(expected) {
					area_expected += area_of_a_px;
				}
				if(actual != expected) {
					error_area += area_of_a_px;
				}
			}
		}
//		std::cout << 10000.0f*error_area << " - " << 10000.0f*expected_area << " - " << 10000.0f*real_expected_area << std::endl;
		coverage_error = (area_expected == 0) ? 0 : error_area / area_expected;
		area_expected_global = M_PI * opt.base_radius * opt.base_radius;
	}

	// we delay the multiplication with two, so that we do not need to divide
	// by two in the computation of the area error
	thickness_plane *= 2.0f;

}

//------------------------------------------------------------------------------

Clustering::Clustering()
{

}

std::vector<Seed> Clustering::getClusterCentersAsSeeds() const
{
	std::vector<Seed> seeds(cluster.size());
	for(std::size_t i=0; i<cluster.size(); i++) {
		const Cluster& c = cluster[i];
		seeds[i] = Seed{c.center.spatial_x(), c.center.spatial_y(), c.center.image_super_radius};
	}
	return seeds;
}

void Clustering::CreatePoints(const slimage::Image3ub& image, const slimage::Image1ui16& depth, const slimage::Image3f& normals)
{
	slimage::Image3f colf(image.width(), image.height());
	// convert to desired color space
	switch(opt.color_space) {
	case ColorSpaces::RGB: {
		slimage::ParallelProcess(image, colf, [](const slimage::It3ub& cub, const slimage::It3f& cf) {
			float r = float(cub[0]) / 255.0f;
			float g = float(cub[1]) / 255.0f;
			float b = float(cub[2]) / 255.0f;
			cf[0] = r;
			cf[1] = g;
			cf[2] = b;
		}, slimage::ThreadingOptions::Single());
	} break;
	case ColorSpaces::HSV: {
		slimage::ParallelProcess(image, colf, [](const slimage::It3ub& cub, const slimage::It3f& cf) {
			float r = float(cub[0]) / 255.0f;
			float g = float(cub[1]) / 255.0f;
			float b = float(cub[2]) / 255.0f;
			Danvil::color_rgb_to_lab(r, g, b, r, g, b);
			r /= 1000.0f;
			g /= 100.0f;
			b /= 100.0f;
			cf[0] = r;
			cf[1] = g;
			cf[2] = b;
		}, slimage::ThreadingOptions::Single());
	} break;
	case ColorSpaces::HN: {
		slimage::ParallelProcess(image, colf, [](const slimage::It3ub& cub, const slimage::It3f& cf) {
			float r = float(cub[0]) / 255.0f;
			float g = float(cub[1]) / 255.0f;
			float b = float(cub[2]) / 255.0f;
			float a = r + g + b;
			if(a > 0.05f) {
				r /= a;
				g /= a;
				b = a * 0.1;
			}
			else {
				r = 0;
				g = 0;
				b = 0;
			}
			cf[0] = r;
			cf[1] = g;
			cf[2] = b;
		}, slimage::ThreadingOptions::Single());
	} break;
	};
	// compute points
	CreatePoints(colf, depth, normals);
}

void Clustering::CreatePoints(const slimage::Image3f& image, const slimage::Image1ui16& depth, const slimage::Image3f& normals)
{
	assert(normals.isNull());

	const float cTempBaseRadius = 0.025f;

	// compute base radius from desired super pixel count
	if(opt.count > 0) {
		opt.base_radius = cTempBaseRadius;
	}

	unsigned int width = image.width();
	unsigned int height = image.height();
	assert(width == depth.width() && height == depth.height());
	if(!normals.isNull()) {
		assert(width == normals.width() && height == normals.height());
		// FIXME this case is disabled!
	}

	points = ImagePoints(width, height);

	slimage::It3f p_col = image.begin();
	slimage::It1ui16 p_depth = depth.begin();
	//const float* p_normals = normals.isNull() ? 0 : normals.begin();

	for(unsigned int y=0; y<height; y++) {
		for(unsigned int x=0; x<width; x++, ++p_col, ++p_depth) {
			Point& p = points(x, y);
			p.color[0] = p_col[0];
			p.color[1] = p_col[1];
			p.color[2] = p_col[2];
			p.depth_i16 = *p_depth;
			if(p.depth_i16 == 0) {
				p.world = Eigen::Vector3f::Zero();
				p.image_super_radius = 0.0f;
				p.gradient = Eigen::Vector2f::Zero();
				p.circularity = 1.0f;
			}
			else {
				float scala = opt.camera.scala(p.depth_i16);
				p.image_super_radius = opt.base_radius * scala;
	//			p.world = opt.camera.unproject(x, y, p.depth_i16);
	//			p.image_super_radius = opt.computePixelScala(p.depth_i16);
	//			p.gradient = LocalDepthGradient(depth, x, y, opt.base_radius, opt.camera);
				p.world = opt.camera.unprojectUsingScala(x, y, scala);
				// FIXME in count mode the gradient is computed using a default radius of 0.02
				// FIXME regardless of the radius chosen later
				p.gradient = LocalDepthGradientUsingScala(depth, x, y, scala, p.image_super_radius, opt.camera);
				p.circularity = Danvil::MoreMath::FastInverseSqrt(p.gradient.squaredNorm() + 1.0f);
				// limit minimal circularity such that the maximum angle is 80 deg
				p.circularity = std::max(p.circularity, 0.173913043f);
			}
		}
	}

	DANVIL_BENCHMARK_START(density)

	// compute desired density
	density = ComputeDepthDensity(points, opt);
	DANVIL_BENCHMARK_STOP(density)

	if(opt.count > 0) {
		// sum density image
		float density_sum = std::accumulate(density.begin(), density.end(), 0.0f, [](float a, float v) { return a + v; });
		float p = static_cast<float>(opt.count) / density_sum;
		std::for_each(density.begin(), density.end(), [p](const slimage::PixelAccess<slimage::Traits<float,1>>& v) { v *= p; });
		// correct base radius
		opt.base_radius = opt.base_radius / std::sqrt(p);
		// correct image_base_radius
		for(unsigned int i=0; i<points.size(); i++) {
			points[i].image_super_radius *= opt.base_radius / cTempBaseRadius;
		}
	}

}

//void Clustering::ComputeSuperpixels(const slimage::Image1f& edges)
//{
//	std::vector<Seed> seeds = FindSeeds();
//	if(!edges.isNull()) {
//		ImproveSeeds(seeds, edges);
//	}
//	ComputeSuperpixels(seeds);
//}

void Clustering::ComputeSuperpixels(const std::vector<Seed>& seeds)
{
	CreateClusters(seeds);
//	std::cout << std::endl << " 0: n=" << cluster.size() << std::endl;
//	for(unsigned int i=0; i<cluster.size(); i++) {
//		std::cout << cluster[i].pixel_ids.size() << " ";
//	}
//	std::cout << std::endl;
	for(unsigned int i=0; i<opt.iterations; i++) {
		MoveClusters();
//		std::cout << i+1 << ": n=" << cluster.size() << std::endl;
//		for(unsigned int j=0; j<cluster.size(); j++) {
//			std::cout << cluster[j].pixel_ids.size() << " ";
//		}
	}
//	std::cout << std::endl;
	if(opt.is_conquer_enclaves) {
		//ConquerEnclaves();
		ConquerMiniEnclaves();
	}
	// delete empty superpixels
	PurgeInvalidClusters();
}

namespace SegmentExtraction
{
	struct Point {
		int x, y;
	};

	struct Segment {
		int label;
		std::vector<unsigned int> points_indices;
		std::set<int> border_labels;
		bool has_nan_neighbor;
	};

	Segment ExtractSegment(const slimage::Image1i& labels, const slimage::Image1i& visited, const Point& start) {
		// we need to temporarily store if border pixels have been visited
		slimage::Image1i border_visited(visited.width(), visited.height(), slimage::Pixel1i{0});
		// init segment
		Segment result;
		int base_label = labels(start.x, start.y);
		result.label = base_label;
		result.has_nan_neighbor = false;
		// do a flood fill by managing a list of "open" pixels which need investigation
		std::queue<Point> open;
		open.push(start);
		while(!open.empty()) {
			// get the next pixel which needs investigation
			Point p = open.front();
			open.pop();
			// skip pixels which are out of bounds
			if(!labels.isValidIndex(p.x, p.y)) {
				continue;
			}
			// do not visit pixels twice in this run
			if(visited(p.x, p.y) == base_label || border_visited(p.x, p.y) == 1) {
				continue;
			}
			// check label of current point
			int current_label = labels(p.x, p.y);
			if(current_label == base_label) {
				// point belongs to label -> investigate neighbors
				unsigned int index = static_cast<unsigned int>(p.x) + static_cast<unsigned int>(p.y) * labels.width();
				result.points_indices.push_back(index);
				// add neighbors
				open.push(Point{p.x-1, p.y});
				open.push(Point{p.x+1, p.y});
				open.push(Point{p.x, p.y-1});
				open.push(Point{p.x, p.y+1});
				// mark pixel as visited in this run (visited is re-used in later runs!)
				visited(p.x, p.y) = base_label;
			}
			else {
				// point does not belong to label -> add to neighbors
				if(current_label == -1) {
					result.has_nan_neighbor = true;
				}
				else {
					result.border_labels.insert(current_label);
				}
				// mark pixel temporarily as border-visited
				border_visited(p.x, p.y) = 1;
			}
		}
		return result;
	}

	std::vector<Segment> FindAllSegments(const slimage::Image1i& labels) {
		std::vector<Segment> segments;
		slimage::Image1i visited(labels.width(), labels.height(), slimage::Pixel1i{-1});
		for(unsigned int y=0; y<labels.height(); y++) {
			for(unsigned int x=0; x<labels.width(); x++) {
				// skip invalid pixels
				if(labels(x,y) == -1) {
					continue;
				}
				// skip pixels which are already visited in their run
				if(visited(x,y) != -1) {
					continue;
				}
				// find segment
				Segment result = ExtractSegment(labels, visited, Point{x, y});
				segments.push_back(result);
			}
		}
		return segments;
	}
}

void Clustering::ConquerEnclaves()
{
	// compute labels for every pixel
	std::vector<int> labels_v = ComputePixelLabels();
	slimage::Image1i labels(width(), height(), slimage::Buffer<int>(width()*height(), labels_v.begin().base()));
	// find all segments (i.e. connected regions with the same label)
	std::vector<SegmentExtraction::Segment> segments = SegmentExtraction::FindAllSegments(labels);
	// count number of segments per label
	std::vector<unsigned int> label_segment_count(cluster.size(), 0);
	for(const SegmentExtraction::Segment& x : segments) {
		assert(x.label < cluster.size());
		assert(x.label != -1);
		label_segment_count[x.label] ++;
	}
	// sort segments by size
	std::sort(segments.begin(), segments.end(),
			[](const SegmentExtraction::Segment& a, const SegmentExtraction::Segment& b) {
					return a.points_indices.size() < b.points_indices.size();
			});
	// remove enclaves and exclaves
	for(SegmentExtraction::Segment& x : segments) {
		// do not delete last segment of a label
		if(label_segment_count[x.label] == 1) {
			continue;
		}
		label_segment_count[x.label] --;
		// can not delete if no neighbor (i.e. only nan neighbors)
		if(x.border_labels.size() == 0) {
			continue;
		}
		// find neighbor with largest segment
		// FIXME possibly use a big enough neighbor which matches color and normal best
		unsigned int neighbor_best_size = 0;
		unsigned int neighbor_best_label = 0;
		for(int nlab : x.border_labels) {
			unsigned int nsize = cluster[nlab].pixel_ids.size();
			if(nsize > neighbor_best_size) {
				neighbor_best_size = nsize;
				neighbor_best_label = nlab;
			}
		}
		// assign segment to one of the neighbors
		std::cout << "Conquering enclave with " << x.points_indices.size() << " points" << std::endl;
		// remove pixels from old cluster
		std::vector<unsigned int>& original_pixel_ids = cluster[x.label].pixel_ids;
		std::sort(original_pixel_ids.begin(), original_pixel_ids.end());
		std::sort(x.points_indices.begin(), x.points_indices.end());
		std::vector<unsigned int> new_pixel_ids(original_pixel_ids.size());
		auto it = std::set_difference(original_pixel_ids.begin(), original_pixel_ids.end(),
				x.points_indices.begin(), x.points_indices.end(),
				new_pixel_ids.begin());
		original_pixel_ids = std::vector<unsigned int>(new_pixel_ids.begin(), it);
		// add pixels to new cluster
		cluster[neighbor_best_label].pixel_ids.insert(cluster[neighbor_best_label].pixel_ids.begin(), x.points_indices.begin(), x.points_indices.end());
	}
}

void Clustering::ConquerMiniEnclaves()
{
	// edge neighbors are checked later
	int offset[4] = {
			-1, +1, -static_cast<int>(width()), +static_cast<int>(width())
	};
	// compute pixel labels
	std::vector<int> labels_v = ComputePixelLabels();
	slimage::Image1i labels(width(), height(), slimage::Buffer<int>(width()*height(), labels_v.begin().base()));
	// iterate over all pixels
	for(unsigned int y=1; y+1<labels.height(); y++) {
		for(unsigned int x=1; x+1<labels.width(); x++) {
			// compute pixel index and get label
			unsigned int index = x + y * labels.width();
			int lab = labels(x,y);
			// skip invalid pixels
			if(lab == -1) {
				continue;
			}
			// check neighbors
			int neighbors[4]; // will contain the labels of neighbor pixels
			bool are_all_different = true; // will be true if all neighbors have a different label
			bool are_all_nan = true; // will be true if all neighbors are invalid
			for(unsigned int i=0; i<4; i++) {
				int x = labels[index + offset[i]];
				neighbors[i] = x;
				are_all_different = are_all_different && (x != lab);
				are_all_nan = are_all_nan && (x == -1);
			};
			// relabel pixel if isolated and at least one neighbors is valid
			if(are_all_different && !are_all_nan) {
				// remove from old cluster
				cluster[lab].removePixel(index);
				// find best new cluster (valid and nearest in world coordinates)
				int best_lab = 0;
				float best_dist = 1e9;
				for(unsigned int i=0; i<4; i++) {
					if(neighbors[i] == -1) {
						continue;
					}
					float d2 = (points[index].world - points[index + offset[i]].world).squaredNorm();
					if(d2 < best_dist) {
						best_lab = neighbors[i];
						best_dist = d2;
					}
				}
				// add to new cluster
				assert(best_lab != -1);
				cluster[best_lab].addPixel(index);
				labels(x,y) = best_lab;
			}
		}
	}
}

std::vector<int> Clustering::ComputePixelLabels() const
{
	std::vector<int> labels(points.size(), -1);
	for(unsigned int j=0; j<cluster.size(); j++) {
		for(unsigned int i : cluster[j].pixel_ids) {
			labels[i] = static_cast<int>(j);
		}
	}
	return labels;
}

slimage::Image1i Clustering::ComputeLabels() const
{
	slimage::Image1i img(width(), height(), slimage::Pixel1i{-1});
	for(unsigned int j=0; j<cluster.size(); j++) {
		for(unsigned int i : cluster[j].pixel_ids) {
			img[i] = static_cast<int>(j);
		}
	}
	return img;
}
void Clustering::CreateClusters(const std::vector<Seed>& seeds)
{
	// create clusters
	cluster.clear();
	cluster.reserve(seeds.size());
	for(const Seed& p : seeds) {
		Cluster c;
//		c.center.valid_ = true;
		c.center.pos[0] = float(p.x);
		c.center.pos[1] = float(p.y);
		c.center.image_super_radius = p.scala;
		// assign points (use only half the radius)
		int R = static_cast<int>(std::ceil(c.center.image_super_radius * 0.35f));
		R = std::min(2, R);
		assert(R >= 0 && "CreateClusters: Invalid radius!");
		c.pixel_ids.reserve((2*R + 1)*(2*R + 1));
		unsigned int xmin = std::max<int>(p.x - R, 0);
		unsigned int xmax = std::min<int>(p.x + R, int(points.width()) - 1);
		unsigned int ymin = std::max<int>(p.y - R, 0);
		unsigned int ymax = std::min<int>(p.y + R, int(points.height()) - 1);
		for(unsigned int yi=ymin; yi<=ymax; yi++) {
			for(unsigned int xi=xmin; xi<=xmax; xi++) {
				unsigned int index = points.index(xi, yi);
				if(points(xi, yi).isInvalid()) {
					// omit invalid points
					continue;
				}
				c.pixel_ids.push_back(index);
			}
		}
		// update center
		if(c.hasPoints()) {
			c.UpdateCenter(points, opt);
			cluster.push_back(c);
		}
	}
}

slimage::Image1f Clustering::ComputeEdges()
{
	const unsigned int width = points.width();
	const unsigned int height = points.height();
	slimage::Image1f edges(width, height);

	// compute edges strength
	slimage::It1f p_edge_begin = edges.begin();
	slimage::ParallelProcess(edges, [this,p_edge_begin,width,height](const slimage::It1f& p_edge) {
		int i = p_edge - p_edge_begin;
		int x = i % width;
		int y = i / width;
		float v = 0.0f;
		if(x < 1 || int(width) <= x+1 || y < 1 || int(height) <= y+1) {
			v = 1e6; // dont want to be here
		}
		else {
			const Point& p0 = points(x, y);
			const Point& px1 = points(x-1, y);
			const Point& px2 = points(x+1, y);
			const Point& py1 = points(x, y-1);
			const Point& py2 = points(x, y+1);
			if(p0.isInvalid() || px1.isInvalid() || px2.isInvalid() || py1.isInvalid() || py2.isInvalid()) {
				v = 1e6; // dont want to be here
			}
			else {
				float dx = Distance(px1, px2);
				float dy = Distance(py1, py2);
				v = dx + dy;
			}
		}
		*p_edge = v;
	}, slimage::ThreadingOptions::Single());

	return edges;
}

void Clustering::ImproveSeeds(std::vector<Seed>& seeds, const slimage::Image1f& edges)
{
	const unsigned int width = points.width();
	const unsigned int height = points.height();

	const int dx8[8] = {-1, -1,  0,  1, 1, 1, 0, -1};
	const int dy8[8] = { 0, -1, -1, -1, 0, 1, 1,  1};

	for(Seed& seed : seeds) {
		int sx = seed.x;
		int sy = seed.y;
		int bestid = sy*width + sx;
		for(unsigned int i=0; i<8; i++) {
			int nx = sx + dx8[i];
			int ny = sy + dy8[i];
			if(nx >= 0 && nx < int(width) && ny >= 0 && ny < int(height)) {
				int nind = ny*width + nx;
				if(edges[nind] < edges[bestid]) {
					bestid = nind;
				}
			}
		}
		seed.x = bestid % width;
		seed.y = bestid / width;
	}
}

void Clustering::PurgeInvalidClusters()
{
	auto it = std::remove_if(cluster.begin(), cluster.end(), [](const Cluster& c) { return !c.hasPoints(); });
	cluster.resize(it - cluster.begin());
}

void Clustering::MoveClusters()
{
	std::vector<float> v_dist(points.size(), 1e9);
	std::vector<int> v_label(points.size(), -1);
	// for each cluster check possible points
	for(unsigned int j=0; j<cluster.size(); j++) {
		const Cluster& c = cluster[j];
		int cx = c.center.spatial_x();
		int cy = c.center.spatial_y();
		const int R = int(c.center.image_super_radius * opt.coverage);
		const unsigned int xmin = std::max(0, cx - R);
		const unsigned int xmax = std::min(int(points.width()-1), cx + R);
		const unsigned int ymin = std::max(0, cy - R);
		const unsigned int ymax = std::min(int(points.height()-1), cy + R);
		//unsigned int pnt_index = points.index(xmin,ymin);
		for(unsigned int y=ymin; y<=ymax; y++) {
			for(unsigned int x=xmin; x<=xmax; x++/*, pnt_index++*/) {
				unsigned int pnt_index = points.index(x, y);
				const Point& p = points[pnt_index];
				if(p.isInvalid()) {
					// omit invalid points
					continue;
				}
				float dist = Distance(p, c.center);
				float& v_dist_best = v_dist[pnt_index];
				if(dist < v_dist_best) {
					v_dist_best = dist;
					v_label[pnt_index] = j;
				}
			}
			//pnt_index -= (xmax - xmin + 1);
			//pnt_index += points.width();
		}
	}
	// delete clusters assignments
	for(Cluster& c : cluster) {
		unsigned int n = c.pixel_ids.size();
		c.pixel_ids.clear();
		c.pixel_ids.reserve(n);
	}
	// assign points to clusters
	for(unsigned int i=0; i<points.size(); i++) {
		int label = v_label[i];
		if(label >= 0) {
			cluster[label].pixel_ids.push_back(i);
		}
	}
	// remove invalid clusters
	PurgeInvalidClusters();
	// update remaining (valid) clusters
//#ifdef CLUSTER_UPDATE_MULTITHREADING
//	{
//		slimage::detail::ThreadPoolManager pool(247); // arg!
//		for(unsigned int i=0; i<clusters.size(); i++) {
//			Cluster* c = &(clusters[i]);
//			pool().schedule(boost::bind(&Cluster::UpdateCenter, c, points, opt.camera));
//		}
//	}
//#else
	for(Cluster& c : cluster) {
		c.UpdateCenter(points, opt);
	}
//#endif
}

std::vector<int> ComputeBorderLabels(unsigned int cid, const Clustering& spc, const slimage::Image1i& labels) {
	const int w = static_cast<int>(spc.width());
	const int h = static_cast<int>(spc.height());
	const int d[4] = { -1, +1, -w, +w };
	std::vector<int> border;
	for(unsigned int pid : spc.cluster[cid].pixel_ids) {
		int x = pid % w;
		int y = pid / w;
		if(1 <= x && x+1 <= w && 1 <= y && y+1 <= h) {
			for(int i=0; i<4; i++) {
				int label = labels[pid + d[i]];
				if(label != cid && label != -1) {
					border.push_back(label);
				}
			}
		}
	}
	return border;
}

/** Computes a list of all pixels which are have label cid and have a face neighbor which has label cjd */
std::vector<unsigned int> ComputeBorderPixelsImpl(unsigned int cid, unsigned int cjd, const Clustering& spc, const slimage::Image1i& labels) {
	const int w = static_cast<int>(spc.width());
	const int h = static_cast<int>(spc.height());
	const int d[4] = { -1, +1, -w, +w };
	std::vector<unsigned int> border;
	for(unsigned int pid : spc.cluster[cid].pixel_ids) {
		int x = pid % w;
		int y = pid / w;
		if(1 <= x && x+1 <= w && 1 <= y && y+1 <= h) {
			for(int i=0; i<4; i++) {
				int label = labels[pid + d[i]];
				if(label == static_cast<int>(cjd)) {
					border.push_back(pid);
				}
			}
		}
	}
	return border;
}

std::vector<std::vector<int> > ComputeBorders(const Clustering& spc) {
	slimage::Image1i labels = spc.ComputeLabels();
	std::vector<std::vector<int> > border_pixels(spc.cluster.size());
	for(unsigned int cid=0; cid<spc.cluster.size(); cid++) {
		border_pixels[cid] = ComputeBorderLabels(cid, spc, labels);
	}
	return border_pixels;
}

std::vector<std::vector<unsigned int> > Clustering::ComputeBorderPixels(const graph::Graph& graph) const
{
	slimage::Image1i labels = ComputeLabels();
	std::vector<std::vector<unsigned int> > borders(graph.edges.size());
	for(unsigned int k=0; k<graph.edges.size(); k++) {
		// compute pixels which are at the border between superpixels e.a and e.b
		unsigned int i = graph.edges[k].a;
		unsigned int j = graph.edges[k].b;
		// superpixel i should have less points than superpixel j
		if(cluster[i].pixel_ids.size() > cluster[j].pixel_ids.size()) {
			std::swap(i,j);
		}
		// find border pixels
		borders[k] = ComputeBorderPixelsImpl(i, j, *this, labels);
	}
	return borders;
}

graph::Graph Clustering::CreateNeighborhoodGraph(NeighborGraphSettings settings) const
{
	std::vector<std::vector<int> > borders = ComputeBorders(*this);
	graph::Graph G;
	G.nodes_ = cluster.size();
	const float node_distance_threshold = settings.max_spatial_distance_mult * opt.base_radius;
	for(unsigned int i=0; i<G.nodes_; i++) {
		for(unsigned int j=i+1; j<G.nodes_; j++) {
			if(settings.cut_by_spatial) {
				float d = (cluster[i].center.world - cluster[j].center.world).norm();
				// only test if distance is smaller than threshold
				if(d > node_distance_threshold) {
					continue;
				}
			}
			// test if superpixels have a common border
			unsigned int cnt_j_in_i = std::count(borders[i].begin(), borders[i].end(), j);
//			unsigned int cnt_i_in_j = std::count(borders[j].begin(), borders[j].end(), i);
//			assert(cnt_j_in_i == cnt_i_in_j);
			float p = static_cast<float>(cnt_j_in_i) / static_cast<float>(std::min(borders[i].size(),borders[j].size()));
			if(p < settings.min_border_overlap) {
				continue;
			}
			// compute cost using color and normal
			float cost;
			if(settings.cost_function == NeighborGraphSettings::SpatialNormalColor) {
				cost = DistanceSpatialColorNormal(cluster[i], cluster[j]);
			}
			else if(settings.cost_function == NeighborGraphSettings::NormalColor) {
				cost = DistanceColorNormal(cluster[i], cluster[j]);
			}
			else {
				BOOST_ASSERT(false);
				cost = 0.0f;
			}
			G.edges.push_back({i,j,cost});
		}
	}
	return G;
}

Clustering::Segmentation Clustering::CreateSegmentation(const graph::Graph& Gn) const
{
	Segmentation seg;
	seg.cluster_graph = Gn;
	// graph segmentation
	seg.segmentation_graph = graph::MinimalSpanningCutting(seg.cluster_graph, opt.segment_threshold, &seg.cluster_labels);
	// remap labels to get a continuous interval of labels
	std::set<unsigned int> unique_labels_set(seg.cluster_labels.begin(), seg.cluster_labels.end());
	std::vector<unsigned int> unique_labels(unique_labels_set.begin(), unique_labels_set.end());
	seg.segment_count = unique_labels.size();
	for(unsigned int& x : seg.cluster_labels) {
		auto it = std::find(unique_labels.begin(), unique_labels.end(), x);
		x = it - unique_labels.begin();
	}
	return seg;
}

ClusterGroupInfo Clustering::ComputeClusterGroupInfo(unsigned int n, float max_thick)
{
	ClusterGroupInfo cgi;
	cgi.hist_thickness = Histogram<float>(n, 0, max_thick);
	cgi.hist_circularity = Histogram<float>(n, 0, 1);
	cgi.hist_area_quotient = Histogram<float>(n, 0.5f, 2.0f);
	cgi.hist_coverage_error = Histogram<float>(n, 0, 0.10f);
	for(const Cluster& c : cluster) {
		cgi.hist_thickness.add(c.thickness);
		cgi.hist_circularity.add(c.circularity);
		cgi.hist_area_quotient.add(c.area_quotient);
		cgi.hist_coverage_error.add(c.coverage_error);
//		std::cout << ci.t << " " << ci.b << " " << ci.a << std::endl;
	}
	return cgi;
}

Clustering ComputeSuperpixels(const slimage::Image3ub& color, const slimage::Image1ui16& depth, const Parameters& opt)
{
	Clustering clustering;
	clustering.opt = opt;
	ComputeSuperpixelsIncremental(clustering, color, depth);
	return clustering;
}

void ComputeSuperpixelsIncremental(Clustering& clustering, const slimage::Image3ub& color, const slimage::Image1ui16& depth)
{
	if(clustering.opt.is_repair_depth) {
		DANVIL_BENCHMARK_START(dasp_repair)
		RepairDepth(depth, color);
		DANVIL_BENCHMARK_STOP(dasp_repair)
	}

	if(clustering.opt.is_smooth_depth) {
		DANVIL_BENCHMARK_START(dasp_smooth)
		SmoothDepth(depth, color);
		DANVIL_BENCHMARK_STOP(dasp_smooth)
	}

	// compute normals only if necessary
	slimage::Image3f normals;
//	if(super_params_ext.weight_normal > 0.0f) {
//		DANVIL_BENCHMARK_START(normals)
//		slimage::Image3f kinect_points = dasp::PointsAndNormals::ComputePoints(kinect_depth, slimage::ThreadingOptions::UsePool(thread_pool_index_));
//		slimage::Image3f kinect_normals = dasp::PointsAndNormals::ComputeNormals(kinect_depth, kinect_points, slimage::ThreadingOptions::UsePool(thread_pool_index_));
//	//	dasp::PointsAndNormals::ComputeNormalsFast(kinect_depth, kinect_points, kinect_normals);
//		DANVIL_BENCHMARK_STOP(normals)
//	}

	// prepare super pixel points
	DANVIL_BENCHMARK_START(dasp_points)
	ImagePoints old_points = clustering.points;
	clustering.CreatePoints(color, depth, normals);
	DANVIL_BENCHMARK_STOP(dasp_points)

	// compute super pixel seeds
	DANVIL_BENCHMARK_START(dasp_seeds)
	clustering.seeds = clustering.FindSeeds(old_points);
//	std::cout << "Seeds: " << seeds.size() << std::endl;
	DANVIL_BENCHMARK_STOP(dasp_seeds)

	// compute super pixel point edges and improve seeds with it
	if(clustering.opt.is_repair_depth) {
		DANVIL_BENCHMARK_START(dasp_improve)
		slimage::Image1f edges = clustering.ComputeEdges();
		clustering.ImproveSeeds(clustering.seeds, edges);
		DANVIL_BENCHMARK_STOP(dasp_improve)
	}

	// compute clusters
	DANVIL_BENCHMARK_START(dasp_clusters)
	clustering.ComputeSuperpixels(clustering.seeds);
	DANVIL_BENCHMARK_STOP(dasp_clusters)
}

//------------------------------------------------------------------------------
}
//------------------------------------------------------------------------------
