
#include "PointDensity.hpp"
#include <pds/Fattal.hpp>

namespace density
{
	template<typename T, typename Fx, typename Fy>
	Eigen::MatrixXf PointDensityImpl(const std::vector<T>& seeds, const Eigen::MatrixXf& target, Fx fx, Fy fy)
	{
		// radius of box in which to average cluster density
		constexpr int RHO_R = 3;
		// range of kernel s.t. 99.9% of mass is covered
		constexpr float cRange = 3.0f*1.10794f;
		constexpr float cMagicSoftener = 0.5f; // 0.62f;
		Eigen::MatrixXf density = Eigen::MatrixXf::Zero(target.rows(), target.cols());
		for(const T& s : seeds) {
			const int sx = std::round(fx(s));
			const int sy = std::round(fy(s));
			// compute point density as average over a box
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
			const float rho = rho_sum / static_cast<float>(rho_num);
			// seed corresponds to a kernel at position (x,y)
			// with sigma = 1/sqrt(pi*rho)
			// i.e. 1/sigma^2 = pi*rho
			// factor pi is already compensated in kernel
			const float sxf = fx(s);
			const float syf = fy(s);
			const float rho_soft = cMagicSoftener * rho;
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
					float delta = rho_soft * pds::fattal::KernelFunctorSquare(rho_soft*d2);
					density(xi, yi) += delta;
				}
			}
		}
		return density;
	}

	Eigen::MatrixXf PointDensity(const std::vector<Eigen::Vector2f>& seeds, const Eigen::MatrixXf& target)
	{
		return PointDensityImpl(seeds, target,
				[](const Eigen::Vector2f& s) { return s[0]; },
				[](const Eigen::Vector2f& s) { return s[1]; }
		);
	}

}
