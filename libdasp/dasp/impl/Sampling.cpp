/*
 * Sampling.cpp
 *
 *  Created on: Mar 25, 2012
 *      Author: david
 */

//------------------------------------------------------------------------------
#include "Sampling.hpp"
#include "BlueNoise.hpp"
#include "Mipmaps.hpp"
#include <Slimage/Parallel.h>
#include <boost/random.hpp>
#include <cmath>
//------------------------------------------------------------------------------
namespace dasp {
//------------------------------------------------------------------------------

slimage::Image1f ComputeDepthDensity(const ImagePoints& points, const Parameters& opt)
{
	slimage::Image1f density(points.width(), points.height());
	for(unsigned int i=0; i<points.size(); i++) {
		const Point& p = points[i];
		/** Estimated number of super pixels at this point
		 * We assume circular superpixels. So the area A of a superpixel at
		 * point location is R*R*pi and the superpixel density is 1/A.
		 * If the depth information is invalid, the density is 0.
		 */
		float cnt = p.isInvalid() ? 0.0f : 1.0f / (M_PI * p.image_super_radius * p.image_super_radius);
		// Additionally the local gradient has to be considered.
		if(opt.gradient_adaptive_density) {
			cnt /= p.circularity;
		}
		density[i] = cnt;
	}
	return density;
}

template<typename T, typename Fxi, typename Fyi, typename Fx, typename Fy>
slimage::Image1f ComputeDepthDensityFromSeeds_impl(const std::vector<T>& seeds, const slimage::Image1f& target, Fxi fxi, Fyi fyi, Fx fx, Fy fy)
{
	// range R of kernel is s.t. phi(x) >= 0.01 * phi(0) for all x <= R
	const float cRange = 1.21f; // BlueNoise::KernelFunctorInverse(0.01f);
	slimage::Image1f density(target.width(), target.height());
	density.fill(slimage::Pixel1f{0.0f});
	for(const T& s : seeds) {
		int sx = fxi(s);
		int sy = fyi(s);
		float sxf = fx(s);
		float syf = fy(s);
		// seed corresponds to a kernel at position (x,y) with sigma = rho(x,y)^(-1/2)
		float rho = target(sx, sy);
//		if(s.x + 1 < int(target.width()) && s.y + 1 < int(target.height())) {
//			rho += target(s.x + 1, s.y) + target(s.x, s.y + 1) + target(s.x + 1, s.y + 1);
//			rho *= 0.25f;
//		}
		// dimension is 2!
		float norm = rho;
		float sigma = 1.0f / std::sqrt(norm);
		float sigma_inv_2 = norm;
		// kernel influence range
		int R = static_cast<int>(std::ceil(cRange * sigma));
		int xmin = std::max<int>(sx - R, 0);
		int xmax = std::min<int>(sx + R, int(target.width()) - 1);
		int ymin = std::max<int>(sy - R, 0);
		int ymax = std::min<int>(sy + R, int(target.height()) - 1);
		for(int yi=ymin; yi<=ymax; yi++) {
			for(int xi=xmin; xi<=xmax; xi++) {
				float d = Square(static_cast<float>(xi) - sxf) + Square(static_cast<float>(yi) - syf);
				float delta = norm * BlueNoise::KernelFunctorSquare(d*sigma_inv_2);
				density(xi, yi) += delta;
			}
		}
	}
	return density;
}

slimage::Image1f ComputeDepthDensityFromSeeds(const std::vector<Seed>& seeds, const slimage::Image1f& target)
{
	return ComputeDepthDensityFromSeeds_impl(seeds, target,
			[](const Seed& s) { return s.x; },
			[](const Seed& s) { return s.y; },
			[](const Seed& s) { return static_cast<float>(s.x); },
			[](const Seed& s) { return static_cast<float>(s.y); });
}

slimage::Image1f ComputeDepthDensityFromSeeds(const std::vector<Eigen::Vector2f>& seeds, const slimage::Image1f& target)
{
	return ComputeDepthDensityFromSeeds_impl(seeds, target,
			[](const Eigen::Vector2f& s) { return static_cast<float>(s[0] + 0.5f); }, // round to nearest integer (assuming positive)
			[](const Eigen::Vector2f& s) { return static_cast<float>(s[1] + 0.5f); },
			[](const Eigen::Vector2f& s) { return s[0]; },
			[](const Eigen::Vector2f& s) { return s[1]; });
}

//------------------------------------------------------------------------------

boost::mt19937 cGlobalRndRng;

void SetRandomNumberSeed(unsigned int x)
{
	cGlobalRndRng.seed(x);
}

std::vector<Seed> FindSeedsGrid(const ImagePoints& points, const Parameters& opt)
{
	unsigned int width = points.width();
	unsigned int height = points.height();
	const float d = std::sqrt(float(width*height) / float(opt.count));
	const unsigned int Nx = (unsigned int)std::ceil(float(width) / d);
	const unsigned int Ny = (unsigned int)std::ceil(float(height) / d);
	const unsigned int Dx = (unsigned int)std::floor(float(width) / float(Nx));
	const unsigned int Dy = (unsigned int)std::floor(float(height) / float(Ny));
	const unsigned int Hx = Dx/2;
	const unsigned int Hy = Dy/2;
	const float S = float(std::max(Dx, Dy));

//	// assume that everything has a distance of 1.5 meters
//	const float cAssumedDistance = 1.5f;
//	unsigned int R = opt.camera.focal / cAssumedDistance * opt.base_radius;
//	unsigned int Dx = R;
//	unsigned int Dy = R;
//	unsigned int Hx = Dx/2;
//	unsigned int Hy = Dy/2;
//	unsigned int Nx = points.width() / Dx;
//	unsigned int Ny = points.height() / Dy;

	// space seeds evently
	std::vector<Seed> seeds;
	seeds.reserve(Nx*Ny);
	for(unsigned int iy=0; iy<Ny; iy++) {
		unsigned int y = Hy + Dy * iy;
		for(unsigned int ix=0; ix<Nx; ix++) {
			unsigned int x = Hx + Dx * ix;
			seeds.push_back(Seed::Dynamic(x, y, S));
		}
	}

	return seeds;
}

bool FindValidSeedPoint(const ImagePoints& points, int& sx0, int& sy0, int range)
{
	if(range == 0) {
		return (sx0 < int(points.width()) && sy0 < int(points.height()) && points(sx0,sy0).isValid());
	}
	// add random offset to add noise
	boost::variate_generator<boost::mt19937&, boost::uniform_int<> > delta(
			cGlobalRndRng, boost::uniform_int<>(-int(range), +int(range)));
	unsigned int trials = 0;
	while(trials < 100) {
		int sx = sx0 + delta();
		int sy = sy0 + delta();
		if(sx < int(points.width()) && sy < int(points.height()) && points(sx,sy).isValid()) {
			sx0 = sx;
			sy0 = sy;
			return true;
		}
		trials++;
	}
	return false;
}

void FindSeedsDepthMipmap_Walk(
		const ImagePoints& points,
		std::vector<Seed>& seeds,
		const std::vector<slimage::Image1f>& mipmaps,
		int level, unsigned int x, unsigned int y)
{
	static boost::uniform_real<float> rnd(0.0f, 1.0f);
	static boost::variate_generator<boost::mt19937&, boost::uniform_real<float> > die(cGlobalRndRng, rnd);

	const slimage::Image1f& mm = mipmaps[level];

	float v = 1.03f * mm(x, y); // TODO nice magic constant

	if(v > 1.0f && level > 1) { // do not access mipmap 0!
		// go down
		FindSeedsDepthMipmap_Walk(points, seeds, mipmaps, level - 1, 2*x,     2*y    );
		FindSeedsDepthMipmap_Walk(points, seeds, mipmaps, level - 1, 2*x,     2*y + 1);
		FindSeedsDepthMipmap_Walk(points, seeds, mipmaps, level - 1, 2*x + 1, 2*y    );
		FindSeedsDepthMipmap_Walk(points, seeds, mipmaps, level - 1, 2*x + 1, 2*y + 1);
	}
	else {
		if(die() <= v)
		{
			unsigned int half = (1 << (level - 1));
			// create seed in the middle
			int sx = static_cast<int>((x << level) + half);
			int sy = static_cast<int>((y << level) + half);
			// add random offset to add noise
			if(FindValidSeedPoint(points, sx, sy, half/2)) { // place near center
				seeds.push_back(Seed::Dynamic(sx, sy, points(sx, sy).image_super_radius));
			}
		}
	}
}

std::vector<Seed> FindSeedsDepthMipmap(const ImagePoints& points, const slimage::Image1f& density, const Parameters& opt)
{
	// compute mipmaps
	std::vector<slimage::Image1f> mipmaps = Mipmaps::ComputeMipmaps(density, 1);
	// now create pixel seeds
	std::vector<Seed> seeds;
	FindSeedsDepthMipmap_Walk(points, seeds, mipmaps, mipmaps.size() - 1, 0, 0);
	return seeds;
}

void FindSeedsDepthMipmapFS_Walk(
		const ImagePoints& points,
		std::vector<Seed>& seeds,
		const std::vector<slimage::Image1f>& mipmaps,
		const std::vector<slimage::Image1f>& carry_mipmaps,
		int level, unsigned int x, unsigned int y)
{
	const slimage::Image1f& mm = mipmaps[level];
	const slimage::Image1f& carry_mm = carry_mipmaps[level];

	// compute density by multiplying percentage with parent total
	float v = mm(x, y) + carry_mm(x, y);

	if(level <= 1 || v <= 1.5f) {
		if(v >= 0.5f) {
			// set seed point in the middel of the cell
			unsigned int half = (1 << (level - 1));
			int sx = static_cast<int>((x << level) + half);
			int sy = static_cast<int>((y << level) + half);
			if(FindValidSeedPoint(points, sx, sy, half/4)) { // place near center
				seeds.push_back(Seed::Dynamic(sx, sy, points(sx, sy).image_super_radius));
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
		bool xp1ok = (x+1 < mm.width());
		bool yp1ok = (y+1 < mm.height());
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
		FindSeedsDepthMipmapFS_Walk(points, seeds, mipmaps, carry_mipmaps, level - 1, 2*x,     2*y    );
		FindSeedsDepthMipmapFS_Walk(points, seeds, mipmaps, carry_mipmaps, level - 1, 2*x,     2*y + 1);
		FindSeedsDepthMipmapFS_Walk(points, seeds, mipmaps, carry_mipmaps, level - 1, 2*x + 1, 2*y    );
		FindSeedsDepthMipmapFS_Walk(points, seeds, mipmaps, carry_mipmaps, level - 1, 2*x + 1, 2*y + 1);
	}
}

std::vector<Seed> FindSeedsDepthMipmapFS(const ImagePoints& points, const slimage::Image1f& density, const Parameters& opt)
{
	// compute mipmaps
	std::vector<slimage::Image1f> mipmaps = Mipmaps::ComputeMipmaps(density, 1);
	std::vector<slimage::Image1f> carry_mipmaps(mipmaps.size());
	for(unsigned int i=1; i<mipmaps.size(); i++) { // HACK: skip first as we never use it
		carry_mipmaps[i] = slimage::Image1f(mipmaps[i].width(), mipmaps[i].height(), slimage::Pixel1f{0.0f});
	}
	// now create pixel seeds
	std::vector<Seed> seeds;
	FindSeedsDepthMipmapFS_Walk(points, seeds, mipmaps, carry_mipmaps, mipmaps.size() - 1, 0, 0);
	return seeds;
}

std::vector<Seed> FindSeedsDepthBlue(const ImagePoints& points, const slimage::Image1f& density, const Parameters& opt)
{
	// compute blue noise points
	std::vector<BlueNoise::Point> pnts = BlueNoise::Compute(density);
	// convert to seeds
	std::vector<Seed> seeds;
	seeds.reserve(pnts.size());
	for(unsigned int i=0; i<pnts.size(); i++) {
		int x = std::round(pnts[i].x);
		int y = std::round(pnts[i].y);
		if(0 <= x && x < int(points.width()) && 0 <= y && y < int(points.height())) {
			seeds.push_back(Seed::Dynamic(x, y, points(x, y).image_super_radius));
		}
	}
	return seeds;
}

std::vector<Seed> FindSeedsDepthFloyd(const ImagePoints& points, const slimage::Image1f& density, const Parameters& opt)
{
	std::vector<Seed> seeds;
	for(unsigned int y=0; y<density.height() - 1; y++) {
		density(1,y) += density(0,y);
		for(unsigned int x=1; x<density.width() - 1; x++) {
			float v = density(x,y);
			if(v >= 0.5f) {
				v -= 1.0f;
				seeds.push_back(Seed::Dynamic(x, y, points(x, y).image_super_radius));
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

// Variante von Floyd-Steinberg. Vorteil: Keine Schlangenlinien in dünn besetzten Bereichen.
std::vector<Seed> FindSeedsDepthFloydExpo(const ImagePoints& points, const slimage::Image1f& density, const Parameters& opt)
{
	// Fehler der nächsten 8 Zeilen in Ringpuffer speichern
	slimage::Image1f ringbuffer( 16 + density.width(), 8 );
	ringbuffer.fill( {0.0f} );

	// Eine schnelle Zufallszahl
	unsigned int crc32 = 0xffffffff;

	// Bild abtasten
	std::vector<Seed> seeds;
	for(unsigned int y=0; y < density.height(); y++)
	{
		float *pRingBuf = ringbuffer.pointer( 8, y % 8 );
		unsigned int x = 0;
		while( x < density.width() )
		{
			// Dichte an dieser Koordinate
			float v = density(x,y);

			// Zielwert einschließlich diffundiertem Fehler
			float err = v + pRingBuf[ x ];
			if( err >= 0.5f) {
				err-= 1.0f;
				seeds.push_back(Seed::Dynamic(x, y, points(x, y).image_super_radius));
			}

			// Diffundierten Fehler aus dem Ringpuffer löschen,
			// damit die Speicherstelle bei einem erneuten Durchlauf durch den Ringpuffer leer ist.
			pRingBuf[ x ] = 0.0f;

			// Bei Dichte unter 7% den Fehler über Radius 2 diffundieren.
			// Bei Dichte unter 4% den Fehler über Radius 4 diffundieren.
			// Bei Dichte unter 1% den Fehler über Radius 8 diffundieren.
			unsigned int radius = 1 << std::max( 0, std::min( 3, 3 - static_cast< int >( 33.33 * ( fabs( v ) + 0.01 ) ) ) );

			// Dafür sorgen daß die Fehler aller Punkte innerhalb des Radius auf die
			// gleiche Koordinate diffundieren. Sonst akkumuliert sich der Fehler nie.
			int DiffusionX = radius * ( x / radius );
			int DiffusionY = radius * ( y / radius );

			// Zufällig in die eine oder andere Richtung diffundieren,
			// um Spuren zu verwischen.
			if( ( crc32 ^ radius ) & 1 )
			{
				ringbuffer( 8 + DiffusionX + radius, ( DiffusionY          ) % 8 ) += 7.0f / 16.0f * err;
				ringbuffer( 8 + DiffusionX - radius, ( DiffusionY + radius ) % 8 ) += 3.0f / 16.0f * err;
				ringbuffer( 8 + DiffusionX         , ( DiffusionY + radius ) % 8 ) += 5.0f / 16.0f * err;
				ringbuffer( 8 + DiffusionX + radius, ( DiffusionY + radius ) % 8 ) += 1.0f / 16.0f * err;
				crc32 = ( crc32 >> 1 ) ^ 0xedb88320;	// Zufallszahl aktualisieren
			}
			else
			{
				ringbuffer( 8 + DiffusionX + radius, ( DiffusionY          ) % 8 ) += 2.0f / 16.0f * err;
				ringbuffer( 8 + DiffusionX - radius, ( DiffusionY + radius ) % 8 ) += 6.0f / 16.0f * err;
				ringbuffer( 8 + DiffusionX         , ( DiffusionY + radius ) % 8 ) += 2.0f / 16.0f * err;
				ringbuffer( 8 + DiffusionX + radius, ( DiffusionY + radius ) % 8 ) += 6.0f / 16.0f * err;
				crc32 >>= 1;	// Zufallszahl aktualisieren
			}

			++x;
		} // for x
	} // for y
	return seeds;
}

void FindSeedsDeltaMipmap_Walk(const ImagePoints& points, std::vector<Seed>& seeds, const std::vector<slimage::Image2f>& mipmaps, int level, unsigned int x, unsigned int y)
{
	static boost::uniform_real<float> rnd(0.0f, 1.0f);
	static boost::variate_generator<boost::mt19937&, boost::uniform_real<float> > die(cGlobalRndRng, rnd);

	const slimage::Image2f& mm = mipmaps[level];

	float v_sum = mm(x, y)[0];
	float v_abs = std::abs(v_sum);// mm(x, y)[1];

	if(v_abs > 1.0f && level > 1) { // do not access mipmap 0!
		// go down
		FindSeedsDeltaMipmap_Walk(points, seeds, mipmaps, level - 1, 2*x,     2*y    );
		FindSeedsDeltaMipmap_Walk(points, seeds, mipmaps, level - 1, 2*x,     2*y + 1);
		FindSeedsDeltaMipmap_Walk(points, seeds, mipmaps, level - 1, 2*x + 1, 2*y    );
		FindSeedsDeltaMipmap_Walk(points, seeds, mipmaps, level - 1, 2*x + 1, 2*y + 1);
	}
	else {
		if(die() < v_abs)
		{
			unsigned int half = (1 << (level - 1));
			// create seed in the middle
			int sx = (x << level) + half;
			int sy = (y << level) + half;
			// add random offset to add noise
			boost::variate_generator<boost::mt19937&, boost::uniform_int<> > delta(
					cGlobalRndRng, boost::uniform_int<>(-int(half/2), +int(half/2)));
			sx += delta();
			sy += delta();

			if(v_sum > 0.0f) {
				// create seed in the middle
				if(sx < int(points.width()) && sy < int(points.height())) {
					float scala = points(sx, sy).image_super_radius;
					Seed s = Seed::Dynamic(sx, sy, scala);
//					std::cout << s.x << " " << s.y << " " << scala << std::endl;
					seeds.push_back(s);
#ifdef CREATE_DEBUG_IMAGES
					slimage::Image3ub debug = slimage::Ref<unsigned char, 3>(sDebugImages["seeds_delta"]);
					debug(sx, sy) = slimage::Pixel3ub{{255,0,0}};
#endif
				}
			}
			else {
				// find nearest
				int best_dist = 1000000000;
				std::size_t best_index = 0;
				for(std::size_t i=0; i<seeds.size(); i++) {
					const Seed& s = seeds[i];
					// do not remove fixed seed points
					if(s.is_fixed) {
						continue;
					}
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
#ifdef CREATE_DEBUG_IMAGES
				slimage::Image3ub debug = slimage::Ref<unsigned char, 3>(sDebugImages["seeds_delta"]);
				debug(seeds[best_index].x, seeds[best_index].y) = slimage::Pixel3ub{{0,255,255}};
#endif
				seeds.erase(seeds.begin() + best_index);
			}
		}
	}
}

std::vector<Seed> FindSeedsDelta(const ImagePoints& points, const std::vector<Seed>& old_seeds, const slimage::Image1f& density_delta, bool delete_small_scala_seeds)
{
	// compute mipmaps
	std::vector<slimage::Image2f> mipmaps = Mipmaps::ComputeMipmapsWithAbs(density_delta, 1);
	// we need to add and delete points!
	std::vector<Seed> seeds = old_seeds;
	FindSeedsDeltaMipmap_Walk(points, seeds, mipmaps, mipmaps.size() - 1, 0, 0);
//	std::cout << "Delta seeds: " << int(seeds.size()) - int(old_seeds.size()) << std::endl;
	// give all seed points the correct scala
	for(Seed& s : seeds) {
		s.scala = points(s.x, s.y).image_super_radius;
	}
	// delete seeds with low scala
	if(delete_small_scala_seeds) {
		std::vector<Seed> ok_size_seeds;
		ok_size_seeds.reserve(seeds.size());
		for(Seed& s : seeds) {
			if(s.scala >= 2.0f) {
				ok_size_seeds.push_back(s);
			}
		}
		return ok_size_seeds;
	}
	else {
		return seeds;
	}
}

std::vector<Seed> FindSeedsDelta(const ImagePoints& points, const std::vector<Seed>& old_seeds, const ImagePoints& old_points, const slimage::Image1f& density_new, const Parameters& opt)
{
#ifdef CREATE_DEBUG_IMAGES
	slimage::Image3ub debug(points.width(), points.height(), {{0,0,0}});
	sDebugImages["seeds_delta"] = slimage::Ptr(debug);
#endif

	// compute old density
	slimage::Image1f density_old = ComputeDepthDensityFromSeeds(old_seeds, density_new);
	// difference
	slimage::Image1f density_delta = density_new - density_old;
	// use function
	return FindSeedsDelta(points, old_seeds, density_delta, true);
}

//------------------------------------------------------------------------------
}
//------------------------------------------------------------------------------

