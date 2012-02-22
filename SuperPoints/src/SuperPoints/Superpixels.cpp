/*
 * Superpixel.cpp
 *
 *  Created on: Jan 26, 2012
 *      Author: david
 */

#include "Superpixels.hpp"
#include "Mipmaps.hpp"
#include "BlueNoise.hpp"
#include "PointsAndNormals.hpp"
#include <Slimage/Paint.hpp>
#include <Danvil/Tools/MoreMath.h>
#include <eigen3/Eigen/Eigenvalues>
#include <boost/random.hpp>

namespace dasp {

void Cluster::UpdateCenter(const ImagePoints& points, const Camera& cam)
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
	center.pos = cam.project(center.world);
	center.depth_i16 = cam.depth(center.world);

	center.normal = FitNormal(pixel_ids, [&points,&center](unsigned int i) { return points[i].world - center.world; });
	if(center.normal[2] > 0.0f) {
		center.normal *= -1.0f;
	}
	center.gradient[0] = - center.normal[0] / center.normal[2];
	center.gradient[1] = - center.normal[1] / center.normal[2];

//	center.normal = points(center.pos).normal;
//	center.gradient = points(center.pos).gradient;

	// scala is not changed!
}

ParametersExt ComputeParameters(const Parameters& opt, unsigned int width, unsigned int height)
{
	ParametersExt opt_ext(opt);
	opt_ext.width = width;
	opt_ext.height = height;
//	float d = std::sqrt(float(opt_ext.width*opt_ext.height) / float(opt_ext.cluster_count));
//	opt_ext.cluster_nx = (unsigned int)std::ceil(float(opt_ext.width) / d);
//	opt_ext.cluster_ny = (unsigned int)std::ceil(float(opt_ext.height) / d);
//	opt_ext.cluster_dx = (unsigned int)std::floor(float(opt_ext.width) / float(opt_ext.cluster_nx));
//	opt_ext.cluster_dy = (unsigned int)std::floor(float(opt_ext.height) / float(opt_ext.cluster_ny));
//	opt_ext.cluster_count = opt_ext.cluster_nx * opt_ext.cluster_ny;
//	opt_ext.radius = std::sqrt(float(opt_ext.cluster_dx*opt_ext.cluster_dx + opt_ext.cluster_dy*opt_ext.cluster_dy));
//	opt_ext.spatial_normalizer = 1.0f / opt_ext.radius;
//	opt_ext.weight_spatial_final = opt.weight_spatial * opt_ext.spatial_normalizer;
	return opt_ext;
}

ImagePoints CreatePoints(
		const slimage::Image3ub& image,
		const slimage::Image1ui16& depth,
		const slimage::Image3f& normals,
		const ParametersExt& opt)
{
	slimage::Image3f colf(image.width(), image.height());
	slimage::ParallelProcess(image, colf, [](const unsigned char* cub, float* cf) {
		cf[0] = float(cub[0]) / 255.0f;
		cf[1] = float(cub[1]) / 255.0f;
		cf[2] = float(cub[2]) / 255.0f;
	});
	return CreatePoints(colf, depth, normals, opt);
}

ImagePoints CreatePoints(
		const slimage::Image3f& image,
		const slimage::Image1ui16& depth,
		const slimage::Image3f& normals,
		const ParametersExt& opt)
{
	unsigned int width = image.width();
	unsigned int height = image.height();
	assert(width == depth.width() && height == depth.height());
	if(!normals.isNull()) {
		assert(width == normals.width() && height == normals.height());
	}

	ImagePoints points(width, height);

	const float* p_col = image.begin();
	const uint16_t* p_depth = depth.begin();
	const float* p_normals = normals.isNull() ? 0 : normals.begin();

	for(unsigned int y=0; y<height; y++) {
		for(unsigned int x=0; x<width; x++, p_col+=3, p_depth++) {
			Point& p = points(x, y);
			p.color[0] = p_col[0];
			p.color[1] = p_col[1];
			p.color[2] = p_col[2];
			p.depth_i16 = *p_depth;
			p.world = opt.camera.unproject(x, y, p.depth_i16);
			p.image_super_radius = opt.computePixelScala(p.depth_i16);
			p.gradient = LocalDepthGradient(depth, x, y, opt.base_radius, opt.camera);
			if(p_normals != 0) {
				p.normal[0] = p_normals[0];
				p.normal[1] = p_normals[1];
				p.normal[2] = p_normals[2];
				p_normals += 3;
			}
			else {
				p.normal[0] = p.gradient[0];
				p.normal[1] = p.gradient[1];
				p.normal[2] = -1.0f;
				p.normal.normalize();
			}
		}
	}

	return points;
}

std::vector<Cluster> ComputeSuperpixels(const ImagePoints& points, const slimage::Image1f& edges, const ParametersExt& opt)
{
	std::vector<Seed> seeds = FindSeeds(points, opt);
	if(!edges.isNull()) {
		ImproveSeeds(seeds, points, edges, opt);
	}
	return ComputeSuperpixels(points, seeds, opt);
}

std::vector<Cluster> ComputeSuperpixels(const ImagePoints& points, const std::vector<Seed>& seeds, const ParametersExt& opt)
{
	std::vector<Cluster> clusters = CreateClusters(seeds, points, opt);
	for(unsigned int i=0; i<opt.iterations; i++) {
		MoveClusters(clusters, points, opt);
	}
	return clusters;
}

std::vector<int> ComputePixelLabels(const std::vector<Cluster>& clusters, const ImagePoints& points)
{
	std::vector<int> labels(points.size(), -1);
	for(unsigned int j=0; j<clusters.size(); j++) {
		for(unsigned int i : clusters[j].pixel_ids) {
			labels[i] = int(j);
		}
	}
	return labels;
}

std::vector<Seed> FindSeedsGrid(const ImagePoints& points, const ParametersExt& opt)
{
//	const float d = std::sqrt(float(opt.width*opt.height) / float(opt.cluster_count));
//	const unsigned int Nx = (unsigned int)std::ceil(float(opt.width) / d);
//	const unsigned int Ny = (unsigned int)std::ceil(float(opt.height) / d);
//	const unsigned int Dx = (unsigned int)std::floor(float(opt.width) / float(Nx));
//	const unsigned int Dy = (unsigned int)std::floor(float(opt.height) / float(Ny));
//	const unsigned int Hx = Dx/2;
//	const unsigned int Hy = Dy/2;
//	const float S = float(std::max(Dx, Dy));

	// assume that everything has a distance of 1.5 meters
	const float cAssumedDistance = 1.5f;
	unsigned int R = opt.camera.focal / cAssumedDistance * opt.base_radius;
	unsigned int Dx = R;
	unsigned int Dy = R;
	unsigned int Hx = Dx/2;
	unsigned int Hy = Dy/2;
	unsigned int Nx = opt.width / Dx;
	unsigned int Ny = opt.height / Dy;

	// space seeds evently
	std::vector<Seed> seeds;
	seeds.reserve(Nx*Ny);
	for(unsigned int iy=0; iy<Ny; iy++) {
		unsigned int y = Hy + Dy * iy;
		for(unsigned int ix=0; ix<Nx; ix++) {
			unsigned int x = Hx + Dx * ix;
			Seed p;
			p.x = x;
			p.y = y;
			p.scala = R;
			seeds.push_back(p);
		}
	}

	return seeds;
}

std::vector<Seed> FindSeedsDepthRandom(const ImagePoints& points, const ParametersExt& opt)
{
	assert(false && "FindSeedsRandom: Not implemented!");
//	static boost::mt19937 rng;
//	constexpr float cCameraFocal = 25.0f;
//	// for each pixel compute number of expected clusters
//	std::vector<float> cdf(points.size());
//	for(unsigned int i=0; i<points.size(); i++) {
//		uint16_t zi = *(depth->begin() + i);
//		float z = 0.001f * float(zi);
//		float v = z * z;
//		cdf[i] = (i == 0) ? v : (v + cdf[i-1]);
//	}
//	float sum = cdf[cdf.size() - 1];
//	boost::uniform_real<float> rnd(0.0f, sum);
//	boost::variate_generator<boost::mt19937&, boost::uniform_real<float> > die(rng, rnd);
//	// randomly pick clusters based on probability
//	std::vector<Seed> seeds;
//	seeds.reserve(opt.cluster_count);
//	while(seeds.size() < opt.cluster_count) {
//		Seed s;
//		float rnd = die();
//		auto it = std::lower_bound(cdf.begin(), cdf.end(), rnd);
//		unsigned int index = it - cdf.begin() - 1;
//		uint16_t zi = *(depth->begin() + index);
//		if(zi == 0) {
//			continue;
//		}
//		float z = 0.001f * float(zi);
//		s.x = index % opt.width;
//		s.y = index / opt.width;
//		s.radius = cCameraFocal / z;
//		seeds.push_back(s);
//	}
//	return seeds;
}

void FindSeedsDepthMipmap_Walk(
		const ImagePoints& points,
		std::vector<Seed>& seeds,
		const std::vector<slimage::Image1f>& mipmaps,
		int level, unsigned int x, unsigned int y)
{
	static boost::mt19937 rng;
	static boost::uniform_real<float> rnd(0.0f, 1.0f);
	static boost::variate_generator<boost::mt19937&, boost::uniform_real<float> > die(rng, rnd);

	const slimage::Image1f& mm = mipmaps[level];

	float v = mm(x, y);

	if(v > 1.0f && level > 1) { // do not access mipmap 0!
		// go down
		FindSeedsDepthMipmap_Walk(points, seeds, mipmaps, level - 1, 2*x,     2*y    );
		FindSeedsDepthMipmap_Walk(points, seeds, mipmaps, level - 1, 2*x,     2*y + 1);
		FindSeedsDepthMipmap_Walk(points, seeds, mipmaps, level - 1, 2*x + 1, 2*y    );
		FindSeedsDepthMipmap_Walk(points, seeds, mipmaps, level - 1, 2*x + 1, 2*y + 1);
	}
	else {
		if(die() < v)
		{
			// create seed in the middle
			Seed s;
			unsigned int half = (1 << (level - 1));
			s.x = (x << level) + half;
			s.y = (y << level) + half;
			if(s.x < int(points.width()) && s.y < int(points.height())) {
				s.scala = points(s.x, s.y).image_super_radius;
//				std::cout << s.x << " " << s.y << " " << s.radius << " " << points(s.x, s.y).scala << " " << points(s.x, s.y).depth << std::endl;
//				if(s.scala > 2.0f) {
					seeds.push_back(s);
//				}
			}
		}
	}
}

slimage::Image1f ComputeDepthDensity(const ImagePoints& points, const ParametersExt& opt)
{
	slimage::Image1f num(points.width(), points.height());
	for(unsigned int i=0; i<points.height(); i++) {
		for(unsigned int j=0; j<points.width(); j++) {
			const Point& p = points(j,i);
			/** Estimated number of super pixels at this point
			 * We assume circular superpixels. So the area A of a superpixel at
			 * point location is R*R*pi and the superpixel density is 1/A.
			 * If the depth information is invalid, the density is 0.
			 */
			float cnt = p.isInvalid() ? 0.0f : 1.0f / (M_PI * p.image_super_radius * p.image_super_radius);
			// Additionally the local gradient has to be considered.
			// We assume local affine projection and compute the distortion
			// factor.
			float lambda = 1.0f;
			if(opt.gradient_adaptive_density) {
				lambda = std::sqrt(p.gradient.squaredNorm() + 1.0f);
			}
//			lambda = std::min(20.0f, lambda); // limit
			num(j,i) = lambda * cnt;
		}
	}
	return num;
}

//slimage::Image1d ComputeDepthDensityDouble(const ImagePoints& points)
//{
//	slimage::Image1d num(points.width(), points.height());
//	for(unsigned int i=0; i<points.size(); i++) {
//		num[i] = double(points[i].estimatedCount());
//	}
//	return num;
//}

std::vector<Seed> FindSeedsDepthMipmap(const ImagePoints& points, const ParametersExt& opt)
{
	// compute estimated number of seeds per pixel
	slimage::Image1f num = ComputeDepthDensity(points, opt);
	// compute mipmaps
	std::vector<slimage::Image1f> mipmaps = Mipmaps::ComputeMipmaps(num, 1);
	// now create pixel seeds
	std::vector<Seed> seeds;
	FindSeedsDepthMipmap_Walk(points, seeds, mipmaps, mipmaps.size() - 1, 0, 0);
	return seeds;
}

std::vector<Seed> FindSeedsDepthBlue(const ImagePoints& points, const ParametersExt& opt)
{
	// compute estimated number of seeds per pixel
	slimage::Image1f num = ComputeDepthDensity(points, opt);
	// compute blue noise points
	std::vector<BlueNoise::Point> pnts = BlueNoise::Compute(num);
	// convert to seeds
	std::vector<Seed> seeds;
	seeds.reserve(pnts.size());
	for(unsigned int i=0; i<pnts.size(); i++) {
		Seed s;
		s.x = std::round(pnts[i].x);
		s.y = std::round(pnts[i].y);
		if(0 <= s.x && s.x < int(points.width()) && 0 <= s.y && s.y < int(points.height())) {
			s.scala = points(s.x, s.y).image_super_radius;
			seeds.push_back(s);
		}
	}
	return seeds;
}

std::vector<Seed> FindSeedsDepthFloyd(const ImagePoints& points, const ParametersExt& opt)
{
	// compute estimated number of seeds per pixel
	slimage::Image1f density = ComputeDepthDensity(points, opt);
	std::vector<Seed> seeds;
	for(unsigned int y=0; y<density.height() - 1; y++) {
		density(1,y) += density(0,y);
		for(unsigned int x=1; x<density.width() - 1; x++) {
			float v = density(x,y);
			if(v >= 0.5f) {
				v -= 1.0f;
				Seed s;
				s.x = x;
				s.y = y;
				s.scala = points(s.x, s.y).image_super_radius;
				seeds.push_back(s);
			}
			density(x+1,y  ) += 7.0f / 16.0f * v;
			density(x-1,y+1) += 3.0f / 16.0f * v;
			density(x  ,y+1) += 5.0f / 16.0f * v;
			density(x+1,y+1) += 1.0f / 16.0f * v;
		}
		// carry over
		density(0, y+1) += density(density.width()-1, y);
	}
	return seeds;
}

std::vector<Seed> FindSeeds(const ImagePoints& points, const ParametersExt& opt)
{
	switch(opt.seed_mode) {
	case SeedModes::EquiDistant:
		return FindSeedsGrid(points, opt);
	case SeedModes::DepthShooting:
		return FindSeedsDepthRandom(points, opt);
	case SeedModes::DepthMipmap:
		return FindSeedsDepthMipmap(points, opt);
	case SeedModes::DepthBlueNoise:
		return FindSeedsDepthBlue(points, opt);
	case SeedModes::DepthFloyd:
		return FindSeedsDepthFloyd(points, opt);
	default:
		assert(false && "FindSeeds: Unkown mode!");
	};
}

std::vector<Cluster> CreateClusters(const std::vector<Seed>& seeds, const ImagePoints& points, const ParametersExt& opt)
{
	// create clusters
	std::vector<Cluster> clusters;
	clusters.reserve(seeds.size());
	for(const Seed& p : seeds) {
		Cluster c;
//		c.center.valid_ = true;
		c.center.pos[0] = float(p.x);
		c.center.pos[1] = float(p.y);
		c.center.image_super_radius = p.scala;
		// assign points (use only half the radius)
		int R = static_cast<int>(std::ceil(c.center.image_super_radius * 0.5f));
		assert(R >= 0 && "CreateClusters: Invalid radius!");
		c.pixel_ids.reserve((2*R + 1)*(2*R + 1));
		unsigned int xmin = std::max<int>(p.x - R, 0);
		unsigned int xmax = std::min<int>(p.x + R, int(opt.width) - 1);
		unsigned int ymin = std::max<int>(p.y - R, 0);
		unsigned int ymax = std::min<int>(p.y + R, int(opt.height) - 1);
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
			c.UpdateCenter(points, opt.camera);
			clusters.push_back(c);
		}
	}
	return clusters;
}

void ComputeEdges(const ImagePoints& points, slimage::Image1f& edges, const ParametersExt& opt, slimage::ThreadingOptions threadopt)
{
	const unsigned int width = points.width();
	const unsigned int height = points.height();

	// compute edges strength
	edges.resize(width, height);
	float* p_edge_begin = edges.begin();
	slimage::ParallelProcess(edges, [p_edge_begin,width,height,&points,&opt](float* p_edge) {
		int i = p_edge - p_edge_begin;
		int x = i % width;
		int y = i / width;
		float v = 0.0f;
		if(x-1 < 0 || int(width) <= x+1 || y-1 < 0 || int(height) <= y+1) {
			v = 1e9; // dont want to be here
		}
		else {
			const Point& px1 = points(x-1, y);
			const Point& px2 = points(x+1, y);
			const Point& py1 = points(x, y-1);
			const Point& py2 = points(x, y+1);
//			if(!px1.valid || !px2.valid || !py1.valid || !py2.valid) {
//				v = 1e9; // dont want to be here
//			}
//			else {
				float dx = Distance(px1, px2, opt);
				float dy = Distance(py1, py2, opt);
				v = dx + dy;
//			}
		}
		*p_edge = v;
	}, threadopt);
}

void ImproveSeeds(std::vector<Seed>& seeds, const ImagePoints& points, const slimage::Image1f& edges, const ParametersExt& opt)
{
	const unsigned int width = points.width();
	const unsigned int height = points.height();

	const int dx8[8] = {-1, -1,  0,  1, 1, 1, 0, -1};
	const int dy8[8] = { 0, -1, -1, -1, 0, 1, 1,  1};

	const float* p_edges = edges.begin();

	for(Seed& seed : seeds) {
		int sx = seed.x;
		int sy = seed.y;
		int bestid = sy*width + sx;
		for(unsigned int i=0; i<8; i++) {
			int nx = sx + dx8[i];
			int ny = sy + dy8[i];
			if(nx >= 0 && nx < int(width) && ny >= 0 && ny < int(height)) {
				int nind = ny*width + nx;
				if(p_edges[nind] < p_edges[bestid]) {
					bestid = nind;
				}
			}
		}
		seed.x = bestid % width;
		seed.y = bestid / width;
	}
}

void MoveClusters(std::vector<Cluster>& clusters, const ImagePoints& points, const ParametersExt& opt)
{
	std::vector<float> v_dist(points.size(), 1e9);
	std::vector<int> v_label(points.size(), -1);
	// for each cluster check possible points
	for(unsigned int j=0; j<clusters.size(); j++) {
		const Cluster& c = clusters[j];
		int cx = c.center.spatial_x();
		int cy = c.center.spatial_y();
		int R = int(c.center.image_super_radius * opt.coverage);
		unsigned int xmin = std::max(0, cx - R);
		unsigned int xmax = std::min(int(points.width()), cx + R);
		unsigned int ymin = std::max(0, cy - R);
		unsigned int ymax = std::min(int(points.height()), cy + R);
		for(unsigned int y=ymin; y<ymax; y++) {
			for(unsigned int x=xmin; x<xmax; x++) {
				unsigned int i = points.index(x,y);
				const Point& p = points[i];
				if(p.isInvalid()) {
					// omit invalid points
					continue;
				}
				float dist = Distance(p, c, opt);
				if(dist < v_dist[i]) {
					v_dist[i] = dist;
					v_label[i] = j;
				}
			}
		}
	}
	// delete clusters assignments
	for(Cluster& c : clusters) {
		unsigned int n = c.pixel_ids.size();
		c.pixel_ids.clear();
		c.pixel_ids.reserve(n);
	}
	// assign points to clusters
	for(unsigned int i=0; i<points.size(); i++) {
		int label = v_label[i];
		if(label >= 0) {
			clusters[label].pixel_ids.push_back(i);
		}
	}
	// update cluster centers and remove invalid clusters
	std::vector<Cluster> clusters_valid;
	clusters_valid.reserve(clusters.size());
	for(unsigned int i=0; i<clusters.size(); i++) {
		Cluster c = clusters[i];
		if(c.hasPoints()) {
			c.UpdateCenter(points, opt.camera);
			clusters_valid.push_back(c);
		}
	}
	clusters = clusters_valid;
}

ClusterInfo ComputeClusterInfo(const Cluster& cluster, const ImagePoints& points, const ParametersExt& opt)
{
	Eigen::Matrix3f A = PointCovariance(cluster.pixel_ids, [&points,&cluster](unsigned int i) { return points[i].world - cluster.center.world; });
	Eigen::SelfAdjointEigenSolver<Eigen::Matrix3f> solver;
	solver.computeDirect(A);
	ClusterInfo infos;
	infos.t = std::abs(solver.eigenvalues()[0]);
	infos.b = std::abs(solver.eigenvalues()[1]);
	infos.a = std::abs(solver.eigenvalues()[2]);
	infos.area = M_PI * infos.a * infos.b;
	float u = infos.b/infos.a;
	infos.eccentricity = std::sqrt(1.0f - u*u);
	infos.radius = std::sqrt(infos.a * infos.b);

	{
		float error_area = 0;
		float expected_area = 0;
		int cx = cluster.center.spatial_x();
		int cy = cluster.center.spatial_y();
		int R = int(2.0f * cluster.center.image_super_radius);
		unsigned int xmin = std::max(0, cx - R);
		unsigned int xmax = std::min(int(points.width()), cx + R);
		unsigned int ymin = std::max(0, cy - R);
		unsigned int ymax = std::min(int(points.height()), cy + R);
		for(unsigned int y=ymin; y<ymax; y++) {
			for(unsigned int x=xmin; x<xmax; x++) {
				unsigned int i = points.index(x,y);
				const Point& p = points[i];
				bool expected = (cluster.center.world - p.world).squaredNorm() < opt.base_radius*opt.base_radius;
				bool actual = std::find(cluster.pixel_ids.begin(), cluster.pixel_ids.end(), i) != cluster.pixel_ids.end();
				float size_of_a_px = p.depth() / opt.camera.focal;
				float area_of_a_px = size_of_a_px*size_of_a_px*std::sqrt(p.gradient.squaredNorm() + 1.0f);
				if(expected != actual) {
					error_area += area_of_a_px;
				}
				if(expected) {
					expected_area += area_of_a_px;
				}
			}
		}
		float real_expected_area = M_PI * opt.base_radius * opt.base_radius;
//		std::cout << 10000.0f*error_area << " - " << 10000.0f*expected_area << " - " << 10000.0f*real_expected_area << std::endl;
		//float expected_area = opt.base_radius*opt.base_radius*M_PI;
//		infos.coverage = error_area / expected_area;
		infos.coverage = error_area / real_expected_area;
	}

	return infos;
}

ClusterGroupInfo ComputeClusterGroupInfo(const std::vector<ClusterInfo>& cluster_info)
{
	ClusterGroupInfo cgi;
	cgi.hist_eccentricity = Histogram<float>(50, 0, 1);
	cgi.hist_radius = Histogram<float>(50, 0, 0.10f);
	cgi.hist_thickness = Histogram<float>(50, 0, 0.01f);
	for(const ClusterInfo& ci : cluster_info) {
		cgi.hist_eccentricity.add(ci.eccentricity);
		cgi.hist_radius.add(ci.radius);
		cgi.hist_thickness.add(ci.t);
//		std::cout << ci.t << " " << ci.b << " " << ci.a << std::endl;
	}
	return cgi;
}

void PlotCluster(const Cluster& cluster, const ImagePoints& points, const slimage::Image3ub& img)
{
	unsigned char c_col_r = 255.0f * cluster.center.color[0];
	unsigned char c_col_g = 255.0f * cluster.center.color[1];
	unsigned char c_col_b = 255.0f * cluster.center.color[2];
	PlotCluster(cluster, points, img, slimage::Pixel3ub{{c_col_r, c_col_g, c_col_b}});
}

void PlotCluster(const Cluster& cluster, const ImagePoints& points, const slimage::Image3ub& img, const slimage::Pixel3ub& color)
{
	assert(cluster.hasPoints());
	// plot all pixels belonging to the cluster in the color of the cluster center
	for(unsigned int i : cluster.pixel_ids) {
		const Point& p = points[i];
		img(p.spatial_x(), p.spatial_y()) = color;
	}
	// plot the cluster center (using some kind of inverse color)
	int cx = cluster.center.spatial_x();
	int cy = cluster.center.spatial_y();
	if(0 <= cx && cx < int(points.width()) && 0 <= cy && cy < int(points.height())) {
		unsigned char* col = img.pointer(cx, cy);
		col[0] = 255 - color[0];
		col[1] = 255 - color[1];
		col[2] = 255 - color[2];
	}

}

void PlotCluster(const std::vector<Cluster>& clusters, const ImagePoints& points, const slimage::Image3ub& img)
{
//	std::vector<Danvil::ColorUB> cols = {
//			Danvil::Color::Red, Danvil::Color::Green, Danvil::Color::Blue, Danvil::Color::Yellow, Danvil::Color::Cyan, Danvil::Color::Magenta, Danvil::Color::Black
//	};
	img.fill(0);
	for(const Cluster& x : clusters) {
		PlotCluster(x, points, img);
	}
}

void PlotClusterCross(const Cluster& cluster, const slimage::Image3ub& img, const ParametersExt& opt)
{
	unsigned char c_col_r = 255.0f * cluster.center.color[0];
	unsigned char c_col_g = 255.0f * cluster.center.color[1];
	unsigned char c_col_b = 255.0f * cluster.center.color[2];
	slimage::Pixel3ub color{{c_col_r,c_col_g,c_col_b}};

	int cx = cluster.center.spatial_x();
	int cy = cluster.center.spatial_y();

	Eigen::Vector2f g0 = cluster.center.gradient;
	float g_norm_sq = g0.squaredNorm();
	if(g_norm_sq == 0.0f) {
		g0 = Eigen::Vector2f::Unit(0);
	}
	else {
		g0 /= std::sqrt(g_norm_sq);
	}
	float sp_0 = cluster.center.image_super_radius;
	float sp_small = sp_0 / std::sqrt(g_norm_sq + 1.0f);
	int p1x = static_cast<int>(sp_small * g0[0]);
	int p1y = static_cast<int>(sp_small * g0[1]);
	int p2x = static_cast<int>(- sp_0 * g0[1]);
	int p2y = static_cast<int>(+ sp_0 * g0[0]);
//	std::cout << cluster.center.gradient.transpose() << " --- " << sp_0 << " --- " << sp_small << std::endl;
//	std::cout << p1x << "," << p1y << " --- " << p2x << "," << p2y << std::endl;
//	slimage::PaintLine(img, cx, cy, cx + p1x, cy + p1y, color);
//	slimage::PaintLine(img, cx, cy, cx + p2x, cy + p2y, color);
	slimage::PaintEllipse(img, cx, cy, p1x, p1y, p2x, p2y, color);

	//	Eigen::Vector2f p3 = opt.camera.project(cluster.center.world + 0.05f*cluster.center.normal);
	//	int p3x = std::round(p3[0]);
	//	int p3y = std::round(p3[1]);
	//	slimage::PaintLine(img, cx, cy, p3x, p3y, color);
}

void PlotClustersCross(const std::vector<Cluster>& clusters, const slimage::Image3ub& img, const ParametersExt& opt)
{
	for(const Cluster& x : clusters) {
		PlotClusterCross(x, img, opt);
	}
}

void PlotEdges(const std::vector<int>& labels, const slimage::Image3ub& img, unsigned int edge_w, unsigned char edge_r, unsigned char edge_g, unsigned char edge_b)
{
	const int dx8[8] = {-1, -1,  0,  1, 1, 1, 0, -1};
	const int dy8[8] = { 0, -1, -1, -1, 0, 1, 1,  1};

	unsigned int width = img.width();
	unsigned int height = img.height();

	std::vector<bool> istaken(width*height, false);

	unsigned int i = 0;
	for(unsigned int y=0; y<height; y++) {
		for(unsigned int x=0; x<width; x++, i++) {
			unsigned int np = 0;
			for(unsigned int k=0; k<8; k++) {
				int sx = int(x) + dx8[k];
				int sy = int(y) + dy8[k];
				if(sx >= 0 && sx < int(width) && sy >= 0 && sy < int(height)) {
					int si = sy*width + sx;
					//if( false == istaken[si] )//comment this to obtain internal contours
					{
						if(labels[i] != labels[si]) {
							np++;
						}
					}
				}
			}
			if(np > edge_w) {
				img(x,y) = slimage::Pixel3ub{{edge_r, edge_g, edge_b}};
				istaken[i] = true;
			}
		}
	}
}

template<typename K, unsigned int CC>
void PlotSeedsImpl(const std::vector<Seed>& seeds, const slimage::Image<K,CC>& img, const slimage::Pixel<K,CC>& color, int size)
{
	for(Seed s : seeds) {
		// round position
		int px = s.x;
		int py = s.y;
		if(px < 0 || int(img.width()) <= px || py < 0 || int(img.height()) <= py) {
			continue;
		}
		if(size == 1) {
			img(px, py) = color;

		}
		else if(size == 2) {
			// paint a star
			//    X
			//   X X
			//    X
			if(1 <= px) {
				img(px-1, py) = color;
			}
			if(px + 1 < int(img.width())) {
				img(px+1, py) = color;
			}
			if(1 <= py) {
				img(px, py-1) = color;
			}
			if(py + 1 < int(img.height())) {
				img(px, py+1) = color;
			}
		}
		else {
			// paint a circle
			//    X
			//   X X
			//  X   X
			//   X X
			//    X
			if(1 <= px && 1 <= py) {
				img(px-1, py-1) = color;
			}
			if(1 <= px && py + 1 < int(img.height())) {
				img(px-1, py+1) = color;
			}
			if(px + 1 < int(img.width()) && 1 <= py) {
				img(px+1, py-1) = color;
			}
			if(px + 1 < int(img.width()) && py + 1 < int(img.height())) {
				img(px+1, py+1) = color;
			}
			if(2 <= px) {
				img(px-2, py) = color;
			}
			if(px + 2 < int(img.width())) {
				img(px+2, py) = color;
			}
			if(2 <= py) {
				img(px, py-2) = color;
			}
			if(py + 2 < int(img.height())) {
				img(px, py+2) = color;
			}
		}
	}
}

void PlotSeeds(const std::vector<Seed>& seeds, const slimage::Image1ub& img, unsigned char grey, int size)
{
	PlotSeedsImpl(seeds, img, slimage::Pixel1ub{grey}, size);
}

void PlotSeeds(const std::vector<Seed>& seeds, const slimage::Image3ub& img, const slimage::Pixel3ub& color, int size)
{
	PlotSeedsImpl(seeds, img, color, size);
}



}
