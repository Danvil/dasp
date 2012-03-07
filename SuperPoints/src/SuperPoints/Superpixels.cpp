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
#include "TreeReduction.hpp"
#define DANVIL_ENABLE_BENCHMARK
#include <Danvil/Tools/Benchmark.h>
#include <Danvil/Tools/MoreMath.h>
#include <eigen3/Eigen/Eigenvalues>
#include <boost/random.hpp>

#define CLUSTER_UPDATE_MULTITHREADING

//------------------------------------------------------------------------------
namespace dasp {
//------------------------------------------------------------------------------

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
	if(center.normal[2] == 0.0f) {
		center.normal[2] = -0.01f;
	}
	if(center.normal[2] > 0.0f) {
		center.normal *= -1.0f;
	}
	center.gradient[0] = - center.normal[0] / center.normal[2];
	center.gradient[1] = - center.normal[1] / center.normal[2];

//	center.normal = points(center.pos).normal;
//	center.gradient = points(center.pos).gradient;

	// scala is not changed!
}

void Cluster::ComputeClusterInfo(const ImagePoints& points, const Parameters& opt)
{
	Eigen::Matrix3f A = PointCovariance(pixel_ids, [this,&points](unsigned int i) { return points[i].world - center.world; });
	Eigen::SelfAdjointEigenSolver<Eigen::Matrix3f> solver;
	solver.computeDirect(A);
	t = std::abs(solver.eigenvalues()[0]);
	b = std::abs(solver.eigenvalues()[1]);
	a = std::abs(solver.eigenvalues()[2]);
	area = M_PI * a * b;
	float u = b/a;
	eccentricity = std::sqrt(1.0f - u*u);
	radius = std::sqrt(a * b);

	{
		float error_area = 0;
		float expected_area = 0;
		int cx = center.spatial_x();
		int cy = center.spatial_y();
		int R = int(2.0f * center.image_super_radius);
		unsigned int xmin = std::max(0, cx - R);
		unsigned int xmax = std::min(int(points.width()), cx + R);
		unsigned int ymin = std::max(0, cy - R);
		unsigned int ymax = std::min(int(points.height()), cy + R);
		for(unsigned int y=ymin; y<ymax; y++) {
			for(unsigned int x=xmin; x<xmax; x++) {
				unsigned int i = points.index(x,y);
				const Point& p = points[i];
				bool expected = (center.world - p.world).squaredNorm() < opt.base_radius*opt.base_radius;
				bool actual = std::find(pixel_ids.begin(), pixel_ids.end(), i) != pixel_ids.end();
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
		coverage = error_area / real_expected_area;
	}

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
	slimage::ParallelProcess(image, colf, [](const unsigned char* cub, float* cf) {
		cf[0] = float(cub[0]) / 255.0f;
		cf[1] = float(cub[1]) / 255.0f;
		cf[2] = float(cub[2]) / 255.0f;
	}, slimage::ThreadingOptions::Single());
	CreatePoints(colf, depth, normals);
}

void Clustering::CreatePoints(const slimage::Image3f& image, const slimage::Image1ui16& depth, const slimage::Image3f& normals)
{
	unsigned int width = image.width();
	unsigned int height = image.height();
	assert(width == depth.width() && height == depth.height());
	if(!normals.isNull()) {
		assert(width == normals.width() && height == normals.height());
		// FIXME this case is disabled!
	}

	points = ImagePoints(width, height);

	const float* p_col = image.begin();
	const uint16_t* p_depth = depth.begin();
//	const float* p_normals = normals.isNull() ? 0 : normals.begin();

	for(unsigned int y=0; y<height; y++) {
		for(unsigned int x=0; x<width; x++, p_col+=3, p_depth++) {
			Point& p = points(x, y);
			p.color[0] = p_col[0];
			p.color[1] = p_col[1];
			p.color[2] = p_col[2];
			p.depth_i16 = *p_depth;
			if(p.depth_i16 == 0) {
				p.world = Eigen::Vector3f::Zero();
				p.image_super_radius = 0.0f;
				p.gradient = Eigen::Vector2f::Zero();
				p.normal = Eigen::Vector3f::Zero();
			}
			else {
				float scala = opt.camera.scala(p.depth_i16);
	//			p.world = opt.camera.unproject(x, y, p.depth_i16);
	//			p.image_super_radius = opt.computePixelScala(p.depth_i16);
	//			p.gradient = LocalDepthGradient(depth, x, y, opt.base_radius, opt.camera);
				p.world = opt.camera.unprojectUsingScala(x, y, scala);
				p.image_super_radius = opt.base_radius * scala;
				p.gradient = LocalDepthGradientUsingScala(depth, x, y, scala, p.image_super_radius, opt.camera);
	//			if(p_normals != 0) {
	//				p.normal[0] = p_normals[0];
	//				p.normal[1] = p_normals[1];
	//				p.normal[2] = p_normals[2];
	//				p_normals += 3;
	//			}
	//			else {
					p.normal[0] = p.gradient[0];
					p.normal[1] = p.gradient[1];
					p.normal[2] = -1.0f;
					p.normal *= Danvil::MoreMath::FastInverseSqrt(p.normal.squaredNorm());
					//p.normal.normalize();
	//			}
			}
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
	for(unsigned int i=0; i<opt.iterations; i++) {
		MoveClusters();
	}
}

std::vector<int> Clustering::ComputePixelLabels()
{
	std::vector<int> labels(points.size(), -1);
	for(unsigned int j=0; j<cluster.size(); j++) {
		for(unsigned int i : cluster[j].pixel_ids) {
			labels[i] = int(j);
		}
	}
	return labels;
}

std::vector<Seed> FindSeedsGrid(const ImagePoints& points, const Parameters& opt)
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
	unsigned int Nx = points.width() / Dx;
	unsigned int Ny = points.height() / Dy;

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

std::vector<Seed> FindSeedsDepthRandom(const ImagePoints& points, const Parameters& opt)
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
			Seed s;
			unsigned int half = (1 << (level - 1));
			// create seed in the middle
			s.x = (x << level) + half;
			s.y = (y << level) + half;
			// add random offset to add noise

			boost::variate_generator<boost::mt19937&, boost::uniform_int<> > delta(
					rng, boost::uniform_int<>(-int(half/2), +int(half/2)));
			s.x += delta();
			s.y += delta();
			if(s.x < int(points.width()) && s.y < int(points.height())) {
				s.scala = points(s.x, s.y).image_super_radius;
//				std::cout << s.x << " " << s.y << " " << s.radius << " " << points(s.x, s.y).scala << " " << points(s.x, s.y).depth << std::endl;
				if(s.scala >= 2.0f) {
					seeds.push_back(s);
				}
			}
		}
	}
}

slimage::Image1f ComputeDepthDensity(const ImagePoints& points, const Parameters& opt)
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

slimage::Image1f ComputeDepthDensityFromSeeds(const std::vector<Seed>& seeds, const slimage::Image1f& target, const Parameters& opt)
{
	// range R of kernel is s.t. phi(x) >= 0.01 * phi(0) for all x <= R
	const float cRange = 1.21f; // BlueNoise::KernelFunctorInverse(0.01f);
	slimage::Image1f density(target.width(), target.height());
	density.fill(0);
	for(const Seed& s : seeds) {
		// seed corresponds to a kernel at position (x,y) with sigma = rho(x,y)^(-1/2)
		float rho = target(s.x, s.y);
		if(s.x + 1 < int(target.width()) && s.y + 1 < int(target.height())) {
			rho += target(s.x + 1, s.y) + target(s.x, s.y + 1) + target(s.x + 1, s.y + 1);
			rho *= 0.25f;
		}
		// dimension is 2!
		float norm = rho;
		float sigma = 1.0f / std::sqrt(norm);
		float sigma_inv_2 = rho;
		// kernel influence range
		int R = static_cast<int>(std::ceil(cRange * sigma));
		int xmin = std::max<int>(s.x - R, 0);
		int xmax = std::min<int>(s.x + R, int(target.width()) - 1);
		int ymin = std::max<int>(s.y - R, 0);
		int ymax = std::min<int>(s.y + R, int(target.height()) - 1);
		for(int yi=ymin; yi<=ymax; yi++) {
			for(int xi=xmin; xi<=xmax; xi++) {
				float d = static_cast<float>(Square(xi - s.x) + Square(yi - s.y));
				float delta = norm * BlueNoise::KernelFunctorSquare(d*sigma_inv_2);
				density(xi, yi) += delta;
			}
		}
	}
	return density;
}

//slimage::Image1d ComputeDepthDensityDouble(const ImagePoints& points)
//{
//	slimage::Image1d num(points.width(), points.height());
//	for(unsigned int i=0; i<points.size(); i++) {
//		num[i] = double(points[i].estimatedCount());
//	}
//	return num;
//}

std::vector<Seed> FindSeedsDepthMipmap(const ImagePoints& points, const Parameters& opt)
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

std::vector<Seed> FindSeedsDepthBlue(const ImagePoints& points, const Parameters& opt)
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

std::vector<Seed> FindSeedsDepthFloyd(const ImagePoints& points, const Parameters& opt)
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

void FindSeedsDeltaMipmap_Walk(const ImagePoints& points, std::vector<Seed>& seeds, const std::vector<slimage::Image1f>& mipmaps, int level, unsigned int x, unsigned int y)
{
	static boost::mt19937 rng;
	static boost::uniform_real<float> rnd(0.0f, 1.0f);
	static boost::variate_generator<boost::mt19937&, boost::uniform_real<float> > die(rng, rnd);

	const slimage::Image1f& mm = mipmaps[level];

	float v = mm(x, y);

	if(std::abs(v) > 1.0f && level > 1) { // do not access mipmap 0!
		// go down
		FindSeedsDeltaMipmap_Walk(points, seeds, mipmaps, level - 1, 2*x,     2*y    );
		FindSeedsDeltaMipmap_Walk(points, seeds, mipmaps, level - 1, 2*x,     2*y + 1);
		FindSeedsDeltaMipmap_Walk(points, seeds, mipmaps, level - 1, 2*x + 1, 2*y    );
		FindSeedsDeltaMipmap_Walk(points, seeds, mipmaps, level - 1, 2*x + 1, 2*y + 1);
	}
	else {
		if(die() < std::abs(v))
		{
			unsigned int half = (1 << (level - 1));
			// create seed in the middle
			int sx = (x << level) + half;
			int sy = (y << level) + half;
			// add random offset to add noise
			boost::variate_generator<boost::mt19937&, boost::uniform_int<> > delta(
					rng, boost::uniform_int<>(-int(half/2), +int(half/2)));
			sx += delta();
			sy += delta();

			if(v > 0.0f) {
				// create seed in the middle
				if(sx < int(points.width()) && sy < int(points.height())) {
					float scala = points(sx, sy).image_super_radius;
					Seed s{sx, sy, scala};
//					std::cout << s.x << " " << s.y << " " << scala << std::endl;
					if(s.scala >= 2.0f) {
						seeds.push_back(s);
					}
				}
			}
			else {
				// find nearest
				int best_dist = 1000000000;
				std::size_t best_index = 0;
				for(std::size_t i=0; i<seeds.size(); i++) {
					const Seed& s = seeds[i];
					int dist =  Square(sx - s.x) + Square(sy - s.y);
					if(dist < best_dist) {
						best_dist = dist;
						best_index = i;
					}
				}
//				auto it = std::min_element(seeds.begin(), seeds.end(), [sx, sy](const Seed& a, const Seed& b) {
//					return Square(sx - a.x) + Square(sy - a.y) < Square(sx - b.x) + Square(sy - b.y);
//				});
				// delete nearest seed
//				seeds.erase(it);
				seeds.erase(seeds.begin() + best_index);
			}
		}
	}
}

std::vector<Seed> FindSeedsDelta(const ImagePoints& points, const std::vector<Seed>& old_seeds, const ImagePoints& old_points, const Parameters& opt)
{
	// compute desired density
	slimage::Image1f density_new = ComputeDepthDensity(points, opt);
	// compute old density
	slimage::Image1f density_old = ComputeDepthDensityFromSeeds(old_seeds, density_new, opt);
	// difference
	slimage::Image1f density_delta = density_new - density_old;
	// compute mipmaps
	std::vector<slimage::Image1f> mipmaps = Mipmaps::ComputeMipmaps(density_delta, 1);
	// we need to add and delete points!
	std::vector<Seed> seeds = old_seeds;
	FindSeedsDeltaMipmap_Walk(points, seeds, mipmaps, mipmaps.size() - 1, 0, 0);
	std::cout << "Delta seeds: " << int(seeds.size()) - int(old_seeds.size()) << std::endl;
	// give all seed points the correct scala
	for(Seed& s : seeds) {
		s.scala = points(s.x, s.y).image_super_radius;
	}
	// delete seeds with low scala
	std::vector<Seed> ok_size_seeds;
	ok_size_seeds.reserve(seeds.size());
	for(Seed& s : seeds) {
		if(s.scala >= 2.0f) {
			ok_size_seeds.push_back(s);
		}
	}
	return ok_size_seeds;
}

std::vector<Seed> Clustering::FindSeeds()
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

std::vector<Seed> Clustering::FindSeeds(const std::vector<Seed>& old_seeds, const ImagePoints& old_points)
{
	if(opt.seed_mode == SeedModes::Delta) {
		return FindSeedsDelta(points, old_seeds, old_points, opt);
	}
	else {
		return FindSeeds();
	}
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
		int R = static_cast<int>(std::ceil(c.center.image_super_radius * 0.5f));
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
			c.UpdateCenter(points, opt.camera);
			cluster.push_back(c);
		}
	}
}

void Clustering::ComputeEdges(slimage::Image1f& edges)
{
	const unsigned int width = points.width();
	const unsigned int height = points.height();

	// compute edges strength
	edges.resize(width, height);
	float* p_edge_begin = edges.begin();
	slimage::ParallelProcess(edges, [this,p_edge_begin,width,height,&opt](float* p_edge) {
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
				float dx = Distance(px1, px2);
				float dy = Distance(py1, py2);
				v = dx + dy;
//			}
		}
		*p_edge = v;
	}, threadopt);
}

void Clustering::ImproveSeeds(std::vector<Seed>& seeds, const slimage::Image1f& edges)
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
		unsigned int pnt_index = points.index(xmin,ymin);
		for(unsigned int y=ymin; y<=ymax; y++) {
			for(unsigned int x=xmin; x<=xmax; x++, pnt_index++) {
				const Point& p = points[pnt_index];
				if(p.isInvalid()) {
					// omit invalid points
					continue;
				}
				float dist = Distance(p, c);
				float& v_dist_best = v_dist[pnt_index];
				if(dist < v_dist_best) {
					v_dist_best = dist;
					v_label[pnt_index] = j;
				}
			}
			pnt_index -= (xmax - xmin + 1);
			pnt_index += points.width();
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
	std::vector<Cluster> clusters_valid;
	clusters_valid.reserve(cluster.size());
	for(unsigned int i=0; i<cluster.size(); i++) {
		const Cluster& c = cluster[i];
		if(c.hasPoints()) {
			clusters_valid.push_back(c);
		}
	}
	cluster = clusters_valid;
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
		c.UpdateCenter(points, opt.camera);
	}
//#endif
}

SuperpixelGraph Clustering::CreateNeighborhoodGraph()
{
	SuperpixelGraph G;
	for(const dasp::Cluster& c : cluster) {
		SuperpixelState s;
		s.x = c.center.spatial_x();
		s.y = c.center.spatial_y();
		s.color = c.center.color;
		s.normal = c.center.normal;
		s.position = c.center.world;
		s.scala = c.center.image_super_radius;
		G.nodes_.push_back(s);
	}
	G.createConnections(3.0f * opt.base_radius);
	return G;
}

std::vector<unsigned int> Clustering::CreateSegments(const SuperpixelGraph& Gn, unsigned int* cnt_label)
{
	// create graph
	Romeo::TreeReduction::Graph graph;
	graph.nodes_ = cluster.size();
	for(unsigned int i=0; i<Gn.node_connections_.size(); i++) {
		for(unsigned int j : Gn.node_connections_[i]) {
			float cost = Clustering::DistanceColorNormal(cluster[i], cluster[j]);
			graph.edges.push_back(Romeo::TreeReduction::Edge{i,j,cost});
		}
	}
	// graph segmentation
	std::vector<unsigned int> labels;
	Romeo::TreeReduction::MinimalSpanningCutting(graph, 10.f, &labels);
	// remap labels
	std::set<unsigned int> unique_labels_set(labels.begin(), labels.end());
	std::vector<unsigned int> unique_labels(unique_labels_set.begin(), unique_labels_set.end());
	*cnt_label = unique_labels.size();
	for(unsigned int& x : labels) {
		auto it = std::find(unique_labels.begin(), unique_labels.end(), x);
		x = it - unique_labels.begin();
	}
	return labels;

//	std::vector<unsigned int> labels(nodes_.size());
//	for(unsigned int i=0; i<labels.size(); i++) {
//		labels[i] = i;
//	}
//	*max_label = nodes_.size() - 1;
//	return labels;
}

ClusterGroupInfo Clustering::ComputeClusterGroupInfo()
{
	ClusterGroupInfo cgi;
	cgi.hist_eccentricity = Histogram<float>(50, 0, 1);
	cgi.hist_radius = Histogram<float>(50, 0, 0.10f);
	cgi.hist_thickness = Histogram<float>(50, 0, 0.01f);
	for(const Cluster& c : cluster) {
		cgi.hist_eccentricity.add(c.eccentricity);
		cgi.hist_radius.add(c.radius);
		cgi.hist_thickness.add(c.t);
//		std::cout << ci.t << " " << ci.b << " " << ci.a << std::endl;
	}
	return cgi;
}

//------------------------------------------------------------------------------
}
//------------------------------------------------------------------------------
