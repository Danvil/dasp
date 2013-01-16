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
#include <boost/random.hpp>
#include <cmath>
//------------------------------------------------------------------------------
namespace dasp {
//------------------------------------------------------------------------------

Eigen::MatrixXf ComputeDepthDensity(const ImagePoints& points, const Parameters& opt)
{
	constexpr float NZ_MIN = 0.174f; // = std::sin(80 deg)

	Eigen::MatrixXf density(points.width(), points.height());
	float* p_density = density.data();
	for(unsigned int i=0; i<points.size(); i++) {
		const Point& p = points[i];
		/** Estimated number of super pixels at this point
		 * We assume circular superpixels. So the area A of a superpixel at
		 * point location is R*R*pi and the superpixel density is 1/A.
		 * If the depth information is invalid, the density is 0.
		 */
		float cnt = 0.0f;
		if(p.is_valid) {
			cnt = 1.0f / (M_PI * p.cluster_radius_px * p.cluster_radius_px);
			// Additionally the local gradient has to be considered.
			if(opt.gradient_adaptive_density) {
				cnt /= std::max(NZ_MIN, p.computeCircularity());
			}
		}
		p_density[i] = cnt;
	}
	return density;
}

template<typename T, typename Fxi, typename Fyi, typename Fx, typename Fy>
Eigen::MatrixXf ComputeDepthDensityFromSeeds_impl(const std::vector<T>& seeds, const Eigen::MatrixXf& target, Fxi fxi, Fyi fyi, Fx fx, Fy fy)
{
	// range R of kernel is s.t. phi(x) >= 0.01 * phi(0) for all x <= R
	constexpr float cRange = 1.21f; // BlueNoise::KernelFunctorInverse(0.01f);
	constexpr float cMagicSoftener = 0.62f;
	Eigen::MatrixXf density = Eigen::MatrixXf::Zero(target.rows(), target.cols());
	for(const T& s : seeds) {
		int sx = fxi(s);
		int sy = fyi(s);
		float sxf = fx(s);
		float syf = fy(s);
		// seed corresponds to a kernel at position (x,y) with sigma = rho(x,y)^(-1/2)
		const float rho = cMagicSoftener*target(sx, sy);
		if(rho == 0) {
			continue;
		}
//		if(s.x + 1 < int(target.width()) && s.y + 1 < int(target.height())) {
//			rho += target(s.x + 1, s.y) + target(s.x, s.y + 1) + target(s.x + 1, s.y + 1);
//			rho *= 0.25f;
//		}
		// kernel influence range
		const int R = static_cast<int>(std::ceil(cRange / std::sqrt(rho)));
		const int xmin = std::max<int>(sx - R, 0);
		const int xmax = std::min<int>(sx + R, int(target.rows()) - 1);
		const int ymin = std::max<int>(sy - R, 0);
		const int ymax = std::min<int>(sy + R, int(target.cols()) - 1);
		for(int yi=ymin; yi<=ymax; yi++) {
			for(int xi=xmin; xi<=xmax; xi++) {
				float d2 = Square(static_cast<float>(xi) - sxf) + Square(static_cast<float>(yi) - syf);
				float delta = rho * BlueNoise::KernelFunctorSquare(rho*d2);
				density(xi, yi) += delta;
			}
		}
	}
	return density;
}

Eigen::MatrixXf ComputeDepthDensityFromSeeds(const std::vector<Seed>& seeds, const Eigen::MatrixXf& target)
{
	return ComputeDepthDensityFromSeeds_impl(seeds, target,
			[](const Seed& s) { return s.x; },
			[](const Seed& s) { return s.y; },
			[](const Seed& s) { return static_cast<float>(s.x); },
			[](const Seed& s) { return static_cast<float>(s.y); });
}

Eigen::MatrixXf ComputeDepthDensityFromSeeds(const std::vector<Eigen::Vector2f>& seeds, const Eigen::MatrixXf& target)
{
	return ComputeDepthDensityFromSeeds_impl(seeds, target,
			[](const Eigen::Vector2f& s) { return static_cast<float>(s[0] + 0.5f); }, // round to nearest integer (assuming positive)
			[](const Eigen::Vector2f& s) { return static_cast<float>(s[1] + 0.5f); },
			[](const Eigen::Vector2f& s) { return s[0]; },
			[](const Eigen::Vector2f& s) { return s[1]; });
}

Eigen::MatrixXf ComputeSaliency(const ImagePoints& points, const Parameters& opt)
{
	const int rows = points.rows();
	const int cols = points.cols();
	const float BR2_INV = 1.0f / (opt.base_radius * opt.base_radius);
	Eigen::MatrixXf saliency_col(rows, cols);
//	Eigen::MatrixXf saliency_norm(rows, cols);
	float* p_saliency_col = saliency_col.data();
//	float* p_saliency_norm = saliency_norm.data();
	for(int y=0; y<cols; ++y) {
		for(int x=0; x<rows; ++x, ++p_saliency_col/*,++p_saliency_norm*/) {
			const Point& p = points(x,y);
			if(!p.is_valid) {
				*p_saliency_col = 0.0f;
//				*p_saliency_norm = 0.0f;
				continue;
			}
			const int r = static_cast<int>(p.cluster_radius_px + 0.5f);
			const int x0 = std::max(x - r, 0);
			const int x1 = std::min(x + r, rows - 1);
			const int y0 = std::max(y - r, 0);
			const int y1 = std::min(y + r, cols - 1);
			// compute mean
			const Eigen::Vector3f mean_col = p.color;
			const Eigen::Vector3f mean_pos = p.position;
//			const Eigen::Vector3f mean_normal = p.normal;
			// compute compression error
			float err_col = 0.0f;
			float err_norm = 0.0f;
			float w_total = 0.0f;
			for(int i=y0; i<=y1; i++) {
				for(int j=x0; j<=x1; j++) {
					const Point& q = points(j,i);
					if(!q.is_valid)
						continue;
					float w = 1.0f / (1.0f + (q.position - mean_pos).squaredNorm() * BR2_INV);
					w_total += w;
					err_col += w * (q.color - mean_col).squaredNorm();
//					err_norm += w * (1.0f - q.normal.dot(mean_normal));
				}
			}
			// write
			*p_saliency_col = std::sqrt(err_col / w_total);
//			*p_saliency_norm = err_norm / w_total;
		}
	}
	// normalize
	{
		const float mean = saliency_col.mean();
		const float min = saliency_col.minCoeff();
		const float max = saliency_col.maxCoeff();
//		std::cout << "color: mean=" << mean << ", min=" << min << ", max=" << max << std::endl;
		saliency_col = (saliency_col.array() - mean)/std::max(mean-min, max-mean);
	}
	// {
	// 	const float mean = saliency_norm.mean();
	// 	const float min = saliency_norm.minCoeff();
	// 	const float max = saliency_norm.maxCoeff();
//	// 	std::cout << "normal: mean=" << mean << ", min=" << min << ", max=" << max << std::endl;
	// 	saliency_norm = (saliency_norm.array() - mean)/std::max(mean-min, max-mean);
	// }
	return saliency_col;// + 0.25f * saliency_norm;
}

void AdaptClusterRadiusBySaliency(ImagePoints& points, const Eigen::MatrixXf& saliency, const Parameters& opt)
{
	const float base = std::max(1.0f, opt.weight_depth);
	auto it_p = points.begin();
	auto it_p_end = points.end();
	auto it_s = saliency.data();
	for(; it_p!=it_p_end; ++it_p, ++it_s) {
		it_p->cluster_radius_px *= std::pow(base, -*it_s);
	}
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
		return (
			sx0 < static_cast<int>(points.width())
			&& sy0 < static_cast<int>(points.height())
			&& points(sx0,sy0).is_valid
		);
	}
	// add random offset to add noise
	boost::variate_generator<boost::mt19937&, boost::uniform_int<> > delta(
			cGlobalRndRng, boost::uniform_int<>(-int(range), +int(range)));
	unsigned int trials = 0;
	while(trials < 100) {
		int sx = sx0 + delta();
		int sy = sy0 + delta();
		if(sx < static_cast<int>(points.width())
			&& sy < static_cast<int>(points.height())
			&& points(sx,sy).is_valid
		) {
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
		const std::vector<Eigen::MatrixXf>& mipmaps,
		int level, unsigned int x, unsigned int y)
{
	static boost::uniform_real<float> rnd(0.0f, 1.0f);
	static boost::variate_generator<boost::mt19937&, boost::uniform_real<float> > die(cGlobalRndRng, rnd);

	const Eigen::MatrixXf& mm = mipmaps[level];

	float v = mm(x, y); // TODO nice magic constant

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
			if(FindValidSeedPoint(points, sx, sy, (3*half)/4)) { // place near center
				seeds.push_back(Seed::Dynamic(sx, sy, points(sx, sy).cluster_radius_px));
			}
		}
	}
}

std::vector<Seed> FindSeedsDepthMipmap(const ImagePoints& points, const Eigen::MatrixXf& density)
{
	// compute mipmaps
	std::vector<Eigen::MatrixXf> mipmaps = Mipmaps::ComputeMipmaps(density, 1);
	// now create pixel seeds
	std::vector<Seed> seeds;
	FindSeedsDepthMipmap_Walk(points, seeds, mipmaps, mipmaps.size() - 1, 0, 0);
	return seeds;
}

void WriteMipmap(std::vector<Eigen::MatrixXf>& mipmaps, int level, int x, int y, float d)
{
	// mipmaps[level](x,y) += d;
	for(int k=level,s=1; k>=1; k--,s*=2) {
		Eigen::MatrixXf& mm = mipmaps[k];
		const float ds = d / static_cast<float>(s*s);
		for(int i=0; i<s; i++) {
			for(int j=0; j<s; j++) {
				mm(x+j,y+i) += ds;
			}
		}
	}
}

void FindSeedsDepthMipmapFS_Walk(
		const ImagePoints& points,
		std::vector<Seed>& seeds,
		std::vector<Eigen::MatrixXf>& mipmaps,
		int level, unsigned int x, unsigned int y)
{
	Eigen::MatrixXf& mm = mipmaps[level];

	// compute density by multiplying percentage with parent total
	float v = mm(x, y);

	if(level <= 1 || v <= 1.5f) {
		if(v >= 0.5f) {
			// set seed point in the middel of the cell
			unsigned int half = (1 << (level - 1));
			int sx = static_cast<int>((x << level) + half);
			int sy = static_cast<int>((y << level) + half);
			if(FindValidSeedPoint(points, sx, sy, (3*half)/4)) { // place near center
				seeds.push_back(Seed::Dynamic(sx, sy, points(sx, sy).cluster_radius_px));
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
//		if(q > 0) {
			float scl = v / q;
			if(xp1ok) 		WriteMipmap(mipmaps, level, x+1, y  , 7.0f*scl);
			if(yp1ok) {
				if(xm1ok) 	WriteMipmap(mipmaps, level, x-1, y+1, 3.0f*scl);
							WriteMipmap(mipmaps, level, x  , y+1, 5.0f*scl);
				if(xp1ok) 	WriteMipmap(mipmaps, level, x+1, y+1, 1.0f*scl);
			}
//		}
	}
	else {
		// go down
		FindSeedsDepthMipmapFS_Walk(points, seeds, mipmaps, level - 1, 2*x,     2*y    );
		FindSeedsDepthMipmapFS_Walk(points, seeds, mipmaps, level - 1, 2*x,     2*y + 1);
		FindSeedsDepthMipmapFS_Walk(points, seeds, mipmaps, level - 1, 2*x + 1, 2*y    );
		FindSeedsDepthMipmapFS_Walk(points, seeds, mipmaps, level - 1, 2*x + 1, 2*y + 1);
	}
}

std::vector<Seed> FindSeedsDepthMipmapFS(const ImagePoints& points, const Eigen::MatrixXf& density)
{
	// compute mipmaps
	std::vector<Eigen::MatrixXf> mipmaps = Mipmaps::ComputeMipmaps(density, 1);
	// now create pixel seeds
	std::vector<Seed> seeds;
	FindSeedsDepthMipmapFS_Walk(points, seeds, mipmaps, mipmaps.size() - 1, 0, 0);
	return seeds;
}

std::vector<Seed> FindSeedsDepthBlue(const ImagePoints& points, const Eigen::MatrixXf& density)
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
			seeds.push_back(Seed::Dynamic(x, y, points(x, y).cluster_radius_px));
		}
	}
	return seeds;
}

std::vector<Seed> FindSeedsDepthFloyd(const ImagePoints& points, const Eigen::MatrixXf& density_inp)
{
	Eigen::MatrixXf density = density_inp;
	std::vector<Seed> seeds;
	for(unsigned int y=0; y<density.cols() - 1; y++) {
		density(1,y) += density(0,y);
		for(unsigned int x=1; x<density.rows() - 1; x++) {
			float v = density(x,y);
			if(v >= 0.5f) {
				v -= 1.0f;
				seeds.push_back(Seed::Dynamic(x, y, points(x, y).cluster_radius_px));
			}
			density(x+1,y  ) += 7.0f / 16.0f * v;
			density(x-1,y+1) += 3.0f / 16.0f * v;
			density(x  ,y+1) += 5.0f / 16.0f * v;
			density(x+1,y+1) += 1.0f / 16.0f * v;
		}
		// carry over
		density(0, y+1) += density(density.rows()-1, y);
	}
	return seeds;
}

// Variante von Floyd-Steinberg. Vorteil: Keine Schlangenlinien in dünn besetzten Bereichen.
std::vector<Seed> FindSeedsDepthFloydExpo(const ImagePoints& points, const Eigen::MatrixXf& density)
{
	// Fehler der nächsten 8 Zeilen in Ringpuffer speichern
	Eigen::MatrixXf ringbuffer( 16 + density.rows(), 8 );
	ringbuffer.fill( {0.0f} );

	// Eine schnelle Zufallszahl
	unsigned int crc32 = 0xffffffff;

	// Bild abtasten
	std::vector<Seed> seeds;
	for(unsigned int y=0; y < density.cols(); y++)
	{
		float *pRingBuf = &ringbuffer( 8, y % 8 );
		unsigned int x = 0;
		while( x < density.rows() )
		{
			// Dichte an dieser Koordinate
			const float v = density(x,y);

			// Zielwert einschließlich diffundiertem Fehler
			float err = v + pRingBuf[ x ];
			if( err >= 0.5f) {
				err-= 1.0f;
				seeds.push_back(Seed::Dynamic(x, y, points(x, y).cluster_radius_px));
			}

			// Diffundierten Fehler aus dem Ringpuffer löschen,
			// damit die Speicherstelle bei einem erneuten Durchlauf durch den Ringpuffer leer ist.
			pRingBuf[ x ] = 0.0f;

			// Bei Dichte über  7% den Fehler über Radius 1 diffundieren.
			// Bei Dichte unter 7% den Fehler über Radius 2 diffundieren.
			// Bei Dichte unter 4% den Fehler über Radius 4 diffundieren.
			// Bei Dichte unter 1% den Fehler über Radius 8 diffundieren.
			const unsigned int LogTable[ 7 ] = { 3, 2, 2, 2, 1, 1, 1 };
			const int t = static_cast< int >( 100.0 * fabs( v ) );
			const unsigned int RadiusLog2 = t >= 7 ? 0 : ( t < 1 ? 3 : LogTable[ t ] );
			const unsigned int radius = 1 << RadiusLog2;

			// Dafür sorgen daß die Fehler aller Punkte innerhalb des Radius auf die
			// gleiche Koordinate diffundieren. Sonst akkumuliert sich der Fehler nie.
			// => Ausrichtung auf ein Vielfaches des Radius
			const int DiffusionX = ( x >> RadiusLog2 ) << RadiusLog2;
			const int DiffusionY = ( y >> RadiusLog2 ) << RadiusLog2;

			// Die nächsten Pixel innerhalb des Radius schneller durchlaufen.
			// Annahme: Die Dichte bleibt konstant, sodaß der Radius nicht geändert werden muß.
			// Dann können alle Fehler auf die gleichen Koordinaten diffundieren.
			++x;
			if( v > 0.5f )
			{
				// Überspringen in dichten Bereichen: Seed-Punkte erzeugen
				while( x < DiffusionX + radius )
				{
					// Fehler der übersprungenen Pixel mitnehmen.
					err += density( x, y ) + pRingBuf[ x ] - 1.0f;
					seeds.push_back(Seed::Dynamic(x, y, points(x, y).cluster_radius_px));
					pRingBuf[ x ] = 0;
					++x;
				}
			}
			else
			{
				// Überspringen in spärlichen Gebieten
				while( x < DiffusionX + radius )
				{
					// Fehler der übersprungenen Pixel mitnehmen.
					err += density( x, y ) + pRingBuf[ x ];
					pRingBuf[ x ] = 0;
					++x;
				}
			}

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
					float scala = points(sx, sy).cluster_radius_px;
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
				if(seeds.size() > 0) {
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
}

std::vector<Seed> FindSeedsDelta(const ImagePoints& points, const std::vector<Seed>& old_seeds, const Eigen::MatrixXf& density_delta, bool delete_small_scala_seeds)
{
	// compute mipmaps
	std::vector<slimage::Image2f> mipmaps = Mipmaps::ComputeMipmapsWithAbs(density_delta, 1);
	// we need to add and delete points!
	std::vector<Seed> seeds = old_seeds;
	FindSeedsDeltaMipmap_Walk(points, seeds, mipmaps, mipmaps.size() - 1, 0, 0);
//	std::cout << "Delta seeds: " << int(seeds.size()) - int(old_seeds.size()) << std::endl;
	// give all seed points the correct scala
	for(Seed& s : seeds) {
		s.scala = points(s.x, s.y).cluster_radius_px;
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

std::vector<Seed> FindSeedsDelta(const ImagePoints& points, const std::vector<Seed>& old_seeds, const ImagePoints& old_points, const Eigen::MatrixXf& density_new)
{
#ifdef CREATE_DEBUG_IMAGES
	slimage::Image3ub debug(points.width(), points.height(), {{0,0,0}});
	sDebugImages["seeds_delta"] = slimage::Ptr(debug);
#endif

	// compute old density
	Eigen::MatrixXf density_old = ComputeDepthDensityFromSeeds(old_seeds, density_new);
	// difference
	Eigen::MatrixXf density_delta = density_new - density_old;
	// use function
	return FindSeedsDelta(points, old_seeds, density_delta, true);
}

//------------------------------------------------------------------------------
}
//------------------------------------------------------------------------------

