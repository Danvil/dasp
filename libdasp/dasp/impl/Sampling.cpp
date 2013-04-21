/*
 * Sampling.cpp
 *
 *  Created on: Mar 25, 2012
 *      Author: david
 */

//------------------------------------------------------------------------------
#include <common/color.hpp>
#include "Sampling.hpp"
#include "../Superpixels.hpp"
#include <pds/BlueNoise.hpp>
#include <pds/Mipmaps.hpp>
#include <Slimage/Paint.hpp>
#include <functional>
#include <boost/random.hpp>
#include <boost/math/constants/constants.hpp>
#include <cmath>

#define CREATE_DEBUG_IMAGES

#ifdef CREATE_DEBUG_IMAGES
	#include <fstream>
	#include <boost/format.hpp>
#endif

//------------------------------------------------------------------------------
namespace dasp {
//------------------------------------------------------------------------------

#ifdef CREATE_DEBUG_IMAGES

template<unsigned int Q>
Eigen::MatrixXf CombineMipmaps(const std::vector<Eigen::MatrixXf>& mm)
{
	Eigen::MatrixXf r = Eigen::MatrixXf::Zero(Q*3*mm[0].rows()/2, Q*mm[0].cols());
	r.block(0, 0, Q*mm[0].rows(), Q*mm[0].cols()) = pds::tools::ScaleUp(mm[0], Q);
	unsigned int y = 0;
	for(unsigned int i=1; i<mm.size(); ++i) {
		r.block(Q*mm[0].rows(), y, Q*mm[i].rows(), Q*mm[i].cols()) = pds::tools::ScaleUp(mm[i], Q);
		y += Q*mm[i].cols();
	}
	return r;
}

void DebugShowMatrix(const Eigen::MatrixXf& mat, const std::string& tag)
{
	const float range = 5000.0f / static_cast<float>((640*480)/25);
	sDebugImages[tag] = slimage::Ptr(
			common::MatrixToImage(mat,
	 			std::bind(&common::IntensityColor, std::placeholders::_1, 0.0f, range)));
}

void DebugWriteMatrix(const Eigen::MatrixXf& mat, const std::string& tag)
{
	std::ofstream ofs(tag);
	for(int i=0; i<mat.rows(); i++) {
		for(int j=0; j<mat.cols(); j++) {
			ofs << mat(i,j);
			if(j+1 != mat.cols()) {
				ofs << "\t";
			}
		}
		if(i+1 != mat.rows()) {
			ofs << "\n";
		}
	}
}

template<unsigned int Q>
void DebugMipmap(const std::vector<Eigen::MatrixXf>& mipmaps, const std::string& tag)
{
	// boost::format fmt(tag + "_%2d");
	// for(std::size_t i=0; i<mipmaps.size(); ++i) {
	// 	const float range = 3000.0f / static_cast<float>(mipmaps[i].rows() * mipmaps[i].cols());
	// 	Eigen::MatrixXf scl = pds::tools::ScaleUp(mipmaps[i], ((Q==2) ? 1 : Q)*(1<<i));
	// 	sDebugImages[(fmt % i).str()] = slimage::Ptr(
	// 		common::MatrixToImage(scl,
	// 			std::bind(&common::IntensityColor, std::placeholders::_1, 0.0f, range)));
	// }
	DebugShowMatrix(CombineMipmaps<Q>(mipmaps), tag);
}

template<unsigned int Q>
void DebugMipmapDelta(const std::vector<Eigen::MatrixXf>& mipmaps, const std::string& tag)
{
	// boost::format fmt(tag + "_%2d");
	// for(std::size_t i=0; i<mipmaps.size(); ++i) {
	// 	const float range = 3000.0f / static_cast<float>(mipmaps[i].rows() * mipmaps[i].cols());
	// 	Eigen::MatrixXf scl = pds::tools::ScaleUp(mipmaps[i], ((Q==2) ? 1 : Q)*(1<<i));
	// 	sDebugImages[(fmt % i).str()] = slimage::Ptr(
	// 		common::MatrixToImage(scl,
	// 			std::bind(&common::PlusMinusColor, std::placeholders::_1, range)));
	// }
	const float range = 2500.0f / static_cast<float>((640*480)/25);
	sDebugImages[tag] = slimage::Ptr(
	 		common::MatrixToImage(CombineMipmaps<Q>(mipmaps),
	 			std::bind(&common::PlusMinusColor, std::placeholders::_1, range)));
}

#endif

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
	const int RHO_R = 3;
	// range R of kernel is s.t. phi(x) >= 0.01 * phi(0) for all x <= R
	constexpr float cRange = 1.21f; // pds::fattal::KernelFunctorInverse(0.01f);
	constexpr float cMagicSoftener = 0.62f;
	Eigen::MatrixXf density = Eigen::MatrixXf::Zero(target.rows(), target.cols());
	for(const T& s : seeds) {
		int sx = fxi(s);
		int sy = fyi(s);
		float sxf = fx(s);
		float syf = fy(s);
		// seed corresponds to a kernel at position (x,y) with sigma = rho(x,y)^(-1/2)
		// const float rho = cMagicSoftener * target(sx, sy);
		float rho_sum = 0.0f;
		unsigned int rho_num = 0;
		for(int i=-RHO_R; i<=+RHO_R; ++i) {
			for(int j=-RHO_R; j<=+RHO_R; ++j) {
				const int sxj = sx + j;
				const int syi = sy + i;
				if( 0 <= sxj && sxj < target.rows() &&
					0 <= syi && syi < target.cols())
				{
					rho_sum += target(sxj, syi);
					rho_num ++;
				}
			}
		}
		if(rho_sum == 0.0f || rho_num == 0) {
			continue;
		}
		const float rho = cMagicSoftener * rho_sum / static_cast<float>(rho_num);
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
				float delta = rho * pds::fattal::KernelFunctorSquare(rho*d2);
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
	// FIXME what is base?
	const float base = 0.5f;
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

void WriteMipmap(std::vector<Eigen::MatrixXf>& mipmaps, int level, int x, int y, float d)
{
	// mipmaps[level](x,y) += d;
	for(int k=level,s=1; k>=0; k--,s*=2) {
		Eigen::MatrixXf& mm = mipmaps[k];
		const float ds = d / static_cast<float>(s*s);
		for(int i=0; i<s; i++) {
			for(int j=0; j<s; j++) {
				mm(x+j,y+i) += ds;
			}
		}
	}
}

template<unsigned int Q>
void FindSeedsDepthMipmapFS_Walk(
		const ImagePoints& points,
		std::vector<Seed>& seeds,
		std::vector<Eigen::MatrixXf>& mipmaps,
		int level, unsigned int x, unsigned int y)
{
	Eigen::MatrixXf& mm = mipmaps[level];

	// compute density by multiplying percentage with parent total
	float v = mm(x, y);

	if(level == 0 || v <= 1.5f) {
		if(v >= 0.5f) {
			// set seed point in the middel of the cell
			unsigned int half = (Q*(1 << level)) / 2;
			int sx = static_cast<int>(Q*(x << level) + half);
			int sy = static_cast<int>(Q*(y << level) + half);
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
		FindSeedsDepthMipmapFS_Walk<Q>(points, seeds, mipmaps, level - 1, 2*x,     2*y    );
		FindSeedsDepthMipmapFS_Walk<Q>(points, seeds, mipmaps, level - 1, 2*x,     2*y + 1);
		FindSeedsDepthMipmapFS_Walk<Q>(points, seeds, mipmaps, level - 1, 2*x + 1, 2*y    );
		FindSeedsDepthMipmapFS_Walk<Q>(points, seeds, mipmaps, level - 1, 2*x + 1, 2*y + 1);
	}
}

std::vector<Seed> FindSeedsDepthMipmapFS(const ImagePoints& points, const Eigen::MatrixXf& density)
{
	// compute mipmaps
	std::vector<Eigen::MatrixXf> mipmaps = pds::tools::ComputeMipmaps(density, 1);
// #ifdef CREATE_DEBUG_IMAGES
// 	DebugMipmap<2>(mipmaps, "mmfs");
// #endif
	// now create pixel seeds
	std::vector<Seed> seeds;
	FindSeedsDepthMipmapFS_Walk<2>(points, seeds, mipmaps, mipmaps.size() - 1, 0, 0);
	return seeds;
}

std::vector<Seed> FindSeedsDepthMipmapFS640(const ImagePoints& points, const Eigen::MatrixXf& density)
{
	// compute mipmaps
	std::vector<Eigen::MatrixXf> mipmaps = pds::tools::ComputeMipmaps640x480(density);
#ifdef CREATE_DEBUG_IMAGES
	DebugMipmap<5>(mipmaps, "mmfs640");
#endif
	// now create pixel seeds
	std::vector<Seed> seeds;
	const unsigned int l0 = mipmaps.size() - 1;
	for(unsigned int y=0; y<mipmaps[l0].cols(); ++y) {
		for(unsigned int x=0; x<mipmaps[l0].rows(); x++) {
			FindSeedsDepthMipmapFS_Walk<5>(points, seeds, mipmaps, l0, x, y);
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

template<unsigned int Q>
void FindSeedsDeltaMipmap_Walk(const ImagePoints& points, std::vector<Seed>& seeds,
	const std::vector<Eigen::MatrixXf>& mipmaps_value,
	const std::vector<Eigen::MatrixXf>& mipmaps_delta,
	const std::vector<Eigen::MatrixXf>& mipmaps_delta_abs,
	int level, unsigned int x, unsigned int y)
{
	constexpr float BREAK_SMOOTH = 2.0f;

	static boost::uniform_real<float> rnd(0.0f, 1.0f);
	static boost::variate_generator<boost::mt19937&, boost::uniform_real<float> > die(cGlobalRndRng, rnd);

	const Eigen::MatrixXf& mm_v = mipmaps_value[level];
	const Eigen::MatrixXf& mm_s = mipmaps_delta[level];
	const Eigen::MatrixXf& mm_a = mipmaps_delta_abs[level];

	const float v_val = mm_v(x, y);
	const float v_sum = mm_s(x, y);
	const float v_sum_abs = std::abs(v_sum);
	const float v_abs = mm_a(x, y);

	if(
		level > 0 // do not access mipmap 0!
		&& v_abs > 1.0f
		&& v_val > BREAK_SMOOTH*BREAK_SMOOTH
		//|| (std::abs(v_sum) - v_abs) / (std::abs(v_sum) + v_abs)
	) {
		// go down
		FindSeedsDeltaMipmap_Walk<Q>(points, seeds, mipmaps_value, mipmaps_delta, mipmaps_delta_abs, level - 1, 2*x,     2*y    );
		FindSeedsDeltaMipmap_Walk<Q>(points, seeds, mipmaps_value, mipmaps_delta, mipmaps_delta_abs, level - 1, 2*x,     2*y + 1);
		FindSeedsDeltaMipmap_Walk<Q>(points, seeds, mipmaps_value, mipmaps_delta, mipmaps_delta_abs, level - 1, 2*x + 1, 2*y    );
		FindSeedsDeltaMipmap_Walk<Q>(points, seeds, mipmaps_value, mipmaps_delta, mipmaps_delta_abs, level - 1, 2*x + 1, 2*y + 1);
	}
	else {
		// const float brok = std::abs(v_sum_abs - v_abs) / (v_sum_abs + v_abs);

		// std::cout << level << "," << x << "," << y << ": t=" << v_val << ", r=" << v_sum << ", |r|=" << v_sum_abs << ", a=" << v_abs << "; q=" << brok << std::endl;

		// if(brok > 0.5f) {
		// 	return;
		// }
		if(die() < v_sum_abs)
		{
			const Eigen::MatrixXf& mmd0 = mipmaps_delta[0];

			const unsigned int x0 = std::min<unsigned int>(mmd0.rows(), (x << level));
			const unsigned int y0 = std::min<unsigned int>(mmd0.cols(), (y << level));
			const unsigned int x1 = std::min<unsigned int>(mmd0.rows(), ((x+1) << level));
			const unsigned int y1 = std::min<unsigned int>(mmd0.cols(), ((y+1) << level));

#ifdef CREATE_DEBUG_IMAGES
			slimage::Image3ub debug = slimage::Ref<unsigned char, 3>(sDebugImages["seeds_delta"]);
			auto color = (v_sum > 0.0f ? slimage::Pixel3ub{{255,0,0}} : slimage::Pixel3ub{{0,255,255}});
			slimage::PaintLine(debug, Q*x0, Q*y0, Q*x1, Q*y0, color);
			slimage::PaintLine(debug, Q*x0, Q*y1, Q*x1, Q*y1, color);
			slimage::PaintLine(debug, Q*x0, Q*y0, Q*x0, Q*y1, color);
			slimage::PaintLine(debug, Q*x1, Q*y0, Q*x1, Q*y1, color);
#endif

			if(v_sum > 0.0f) {
				// find coordinate in cell where delta density is minimal
				const auto& b = mmd0.block(x0, y0, 1 << level, 1 << level);
				unsigned int best_j=-1, best_i=-1;
				float best_val = -1000000.0f;
				for(unsigned int i=0; i<b.cols(); ++i) {
					for(unsigned int j=0; j<b.rows(); ++j) {
						float val = b(j,i);
						if(val > best_val) {
							best_j = j;
							best_i = i;
							best_val = val;
						}
					}
				}
				// add seed to middle of cell
				if(best_i != -1 && best_j != -1) {
					int sx = Q*(x0 + best_j) + Q/2;
					int sy = Q*(y0 + best_i) + Q/2;
					float scala = points(sx, sy).cluster_radius_px;
					Seed s = Seed::Dynamic(sx, sy, scala);
//					std::cout << s.x << " " << s.y << " " << scala << std::endl;
					seeds.push_back(s);
#ifdef CREATE_DEBUG_IMAGES
					for(unsigned int i=0; i<2; i++) {
						for(unsigned int j=0; j<2; j++) {
							debug(sx+j, sy+i) = slimage::Pixel3ub{{255,0,0}};
						}
					}
#endif
				}
			}
			else {
				if(seeds.empty()) {
					return;
				}
				// find nearest
				float best_val = +1000000.0f;
				std::size_t best_index = -1;
				for(std::size_t i=0; i<seeds.size(); i++) {
					const Seed& s = seeds[i];
					int sx = s.x / Q;
					int sy = s.y / Q;
					// do not remove fixed seed points
					if(s.is_fixed) {
						continue;
					}
					if(sx < x0 || x1 <= sx || sy < y0 || y1 <= sy) {
						continue;
					}
					float val = mmd0(sx, sy);
					if(val < best_val) {
						best_index = i;
						best_val = val;
					}
				}
				// delete nearest seed
				if(best_index != -1) {
#ifdef CREATE_DEBUG_IMAGES
					unsigned int sx = seeds[best_index].x;
					unsigned int sy = seeds[best_index].y;
					for(unsigned int i=0; i<2; i++) {
						for(unsigned int j=0; j<2; j++) {
							debug(sx+j, sy+i) = slimage::Pixel3ub{{0,255,255}};
						}
					}
#endif
					seeds.erase(seeds.begin() + best_index);
				}
			}
		}
	}
}

std::vector<Seed> FindSeedsDelta(const ImagePoints& points, const std::vector<Seed>& old_seeds, const Eigen::MatrixXf& density_old, const Eigen::MatrixXf& density_new, bool delete_small_scala_seeds)
{
	// difference
	Eigen::MatrixXf density_delta = density_new - density_old;
	// compute mipmaps
	std::vector<Eigen::MatrixXf> mm_v = pds::tools::ComputeMipmaps640x480(density_new);
	std::vector<Eigen::MatrixXf> mm_dv = pds::tools::ComputeMipmaps640x480(density_delta);
	std::vector<Eigen::MatrixXf> mm_da = pds::tools::ComputeMipmaps640x480(density_delta.cwiseAbs());
#ifdef CREATE_DEBUG_IMAGES
	DebugMipmap<5>(mm_v, "mm_v");
	DebugMipmapDelta<5>(mm_dv, "mm_dv");
	DebugMipmap<5>(mm_da, "mm_da");
#endif
	// we need to add and delete points!
	std::vector<Seed> seeds = old_seeds;
	const unsigned int l0 = mm_dv.size() - 1;
	for(unsigned int y=0; y<mm_dv[l0].cols(); ++y) {
		for(unsigned int x=0; x<mm_dv[l0].rows(); x++) {
			FindSeedsDeltaMipmap_Walk<5>(points, seeds, mm_v, mm_dv, mm_da, l0, x, y);
		}
	}
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
	// use function
	return FindSeedsDelta(points, old_seeds, density_old, density_new, true);
}

//------------------------------------------------------------------------------
}
//------------------------------------------------------------------------------

