
#include "Density.hpp"
#include "Fattal.hpp"

namespace pds
{
	template<typename T, typename Fxi, typename Fyi, typename Fx, typename Fy>
	Eigen::MatrixXf ComputeDepthDensityFromSeeds_impl(const std::vector<T>& seeds, const Eigen::MatrixXf& target, Fxi fxi, Fyi fyi, Fx fx, Fy fy)
	{
		const int RHO_R = 3;
		// range R of kernel is s.t. phi(x) >= 0.01 * phi(0) for all x <= R
		constexpr float cRange = 1.21f; // pds::fattal::KernelFunctorInverse(0.01f);
		constexpr float cMagicSoftener = 0.5f;// 0.62f;
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
					float dx = static_cast<float>(xi) - sxf;
					float dy = static_cast<float>(yi) - syf;
					float d2 = dx*dx + dy*dy;
					float delta = rho * pds::fattal::KernelFunctorSquare(rho*d2);
					density(xi, yi) += delta;
				}
			}
		}
		return density;
	}

	Eigen::MatrixXf PointDensity(const std::vector<Eigen::Vector2f>& seeds, const Eigen::MatrixXf& target)
	{
		return ComputeDepthDensityFromSeeds_impl(seeds, target,
				[](const Eigen::Vector2f& s) { return (int)std::round(s[0]); }, // round to nearest integer (assuming positive)
				[](const Eigen::Vector2f& s) { return (int)std::round(s[1]); },
				[](const Eigen::Vector2f& s) { return s[0]; },
				[](const Eigen::Vector2f& s) { return s[1]; });
	}

}
