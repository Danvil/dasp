/*
 * BlueNoise.cpp
 *
 *  Created on: Feb 6, 2012
 *      Author: david
 */

#include "Fattal.hpp"
#include "Mipmaps.hpp"
#include <boost/random.hpp>
//----------------------------------------------------------------------------//
namespace pds {
namespace fattal {
//----------------------------------------------------------------------------//

float EnergyApproximation(const std::vector<Point>& pnts, float x, float y)
{
	float sum = 0.0f;
	for(const Point& p : pnts) {
		float dx = p.x - x;
		float dy = p.y - y;
//			float d = std::sqrt(dx*dx + dy*dy);
//			if(d < KernelRange * p.scale) {
//				float k_val = KernelFunctor(d / p.scale);
		float d2 = dx*dx + dy*dy;
		float scl = p.scale*p.scale;
		if(d2 < KernelRange * KernelRange * scl) {
			float k_val = KernelFunctorSquare(d2 / scl);
			sum += p.weight * ScalePowerD(p.scale) * k_val;
		}
	}
	return sum;
}

float Energy(const std::vector<Point>& pnts, const Eigen::MatrixXf& density)
{
	float error = 0.0f;
	for(unsigned int y=0; y<density.cols(); y++) {
		for(unsigned int x=0; x<density.rows(); x++) {
			float px = float(x);
			float py = float(y);
			float a = EnergyApproximation(pnts, px, py);
			float roh = density(x, y);
			error += std::abs(a - roh);
		}
	}
	return error;
}

float EnergyDerivative(const std::vector<Point>& pnts, const Eigen::MatrixXf& density, unsigned int i, float& result_dE_x, float& result_dE_y)
{
//	result_dE_x = 0.0f;
//	result_dE_y = 0.0f;
//	return;

	float dE_x = 0.0f;
	float dE_y = 0.0f;
	float px = pnts[i].x;
	float py = pnts[i].y;
	float ps = pnts[i].scale;
//		float ps_scl = 1.0f / ps;
	float ps_scl = 1.0f / (ps * ps);
	// find window (points outside the window do not affect the kernel)
	// range = sqrt(ln(1/eps)/pi)
	constexpr float cMaxRange = 1.482837414f; // eps = 0.001
	float radius = cMaxRange * ps;
	float x_min = std::max(0, int(std::floor(px - radius)));
	float x_max = std::min(int(density.rows()) - 1, int(std::ceil(px + radius)));
	float y_min = std::max(0, int(std::floor(py - radius)));
	float y_max = std::min(int(density.cols()) - 1, int(std::ceil(py + radius)));
	// sum over window
	for(unsigned int y=y_min; y<=y_max; y++) {
		for(unsigned int x=x_min; x<=x_max; x++) {
			float ux = float(x);
			float uy = float(y);
			float dx = ux - px;
			float dy = uy - py;
//			float k_arg = std::sqrt(dx*dx + dy*dy) * ps_scl;
//			float k_val = KernelFunctor(k_arg);
			float k_arg_square = (dx*dx + dy*dy) * ps_scl;
			float k_val = KernelFunctorSquare(k_arg_square);
			float apx = EnergyApproximation(pnts, ux, uy);
			float roh = density(x, y);
			if(apx < roh) {
				k_val = -k_val;
			}
			//k_val = 0.0f;
			dE_x += k_val * dx;
			dE_y += k_val * dy;
		}
	}
	float A = 2.0f * cPi / std::pow(ps, float(D + 1));
	result_dE_x = A * dE_x;
	result_dE_y = A * dE_y;
	return radius;
}

std::vector<Point> PlacePoints(const Eigen::MatrixXf& density, unsigned int p)
{
	// access original index in a random order
	std::vector<unsigned int> indices(density.size());
	for(unsigned int i=0; i<indices.size(); i++) {
		indices[i] = i;
	}
	std::random_shuffle(indices.begin(), indices.end());

	/*std::cout << "{" << std::endl;
	for(unsigned int y=0; y<density.height(); y++) {
		std::cout << "{";
		for(unsigned int x=0; x<density.width(); x++) {
			std::cout << density(x,y);
			if(x + 1 < density.width()) {
				std::cout << ", ";
			}
		}
		std::cout << "}";
		if(y + 1 < density.height()) {
			std::cout << ",";
		}
		std::cout << std::endl;
	}
	std::cout << "}" << std::endl;*/

	// compute points
	std::vector<Point> pnts;
	pnts.reserve(indices.size());
	// compute current error in density
	float error_current = Energy(pnts, density);
//	std::cout << "INITIAL ERROR: " << error_current << std::endl;
	// try add kernel points
	for(unsigned int i : indices) {
		float roh = density.data()[i];
		if(roh == 0) {
//				std::cout << i << " roh is 0!" << std::endl;
			continue;
		}
		Point u;
		u.x = float(i % density.rows());
		u.y = float(i / density.rows());
		int q = p - (roh < 1 ? 0 : std::ceil(std::log2(roh) / float(D)));
		u.weight = float(1 << (D*(p-q)));
		u.scale = KernelScaleFunction(roh, u.weight);
		// try to add
		pnts.push_back(u);
		// check if the points reduced the energy
		float error_new = Energy(pnts, density);
		if(error_new > error_current) {
			// reject
			pnts.pop_back();
//			std::cout << u.x << " " << u.y << " " << u.weight << " " << error_new << " REJECTED" << std::endl;
		}
		else {
			error_current = error_new;
//			std::cout << u.x << " " << u.y << " " << u.weight << " " << error_new << std::endl;
		}
	}
	return pnts;
}

void Refine(std::vector<Point>& points, const Eigen::MatrixXf& density, unsigned int iterations)
{
	static boost::mt19937 rng;
	static boost::normal_distribution<float> rnd(0.0f, 1.0f); // standard normal distribution
	static boost::variate_generator<boost::mt19937&, boost::normal_distribution<float> > die(rng, rnd);
	constexpr float dt = 0.2f;
	constexpr float T = 0.2f;
	float r_min = 1e9;
	float r_max = 0;
	for(unsigned int k=0; k<iterations; k++) {
		for(unsigned int i=0; i<points.size(); i++) {
			Point& p = points[i];
			if(p.scale > cMaxRefinementScale) {
				// omit low frequency kernels
				continue;
			}
			// dx = -dt*s/2*dE + sqrt(T*s*dt)*R
			float c0 = dt * p.scale;
			float cA = c0 * 0.5f;
			float dx, dy;
			float R = EnergyDerivative(points, density, i, dx, dy);
			r_min = std::min(R, r_min);
			r_max = std::max(R, r_max);
			p.x -= cA * dx;
			p.y -= cA * dy;
			if(T > 0.0f) {
				float cB = std::sqrt(T * c0);
				p.x += cB * die();
				p.y += cB * die();
			}
		}
	}
//	std::cout << r_min << " " << r_max << std::endl;
}

std::vector<Point> Split(const std::vector<Point>& points, const Eigen::MatrixXf& density, bool& result_added)
{
	std::vector<Point> pnts_new;
	result_added = false;
	for(Point u : points) {
		if(u.weight > 1.0f) {
			result_added = true;
			u.x *= 2.0f;
			u.y *= 2.0f;
			u.weight /= float(1 << D);
			constexpr float A = 0.70710678f;
			constexpr float Delta[4][2] = {
					{-A, -A}, {+A, -A}, {-A, +A}, {+A, +A}
			};
			for(unsigned int i=0; i<4; i++) {
				Point ui = u;
				ui.x += u.scale * Delta[i][0];
				ui.y += u.scale * Delta[i][1];
				float roh = ZeroBorderAccess(density, int(ui.x), int(ui.y));
				if(roh > 0) {
					ui.scale = KernelScaleFunction(roh, ui.weight);
					pnts_new.push_back(ui);
				}
			}
		}
		else {
			u.x *= 2.0f;
			u.y *= 2.0f;
			u.weight = 1.0f;
			float roh = ZeroBorderAccess(density, int(u.x), int(u.y));
			if(roh > 0) {
				u.scale = KernelScaleFunction(roh, u.weight);
				pnts_new.push_back(u);
			}
		}
	}
	return pnts_new;
}

std::vector<Point> Compute(const Eigen::MatrixXf& density, unsigned int max_steps)
{
	// compute mipmaps
	std::vector<Eigen::MatrixXf> mipmaps = tools::ComputeMipmaps(density, 16);
	int p = int(mipmaps.size()) - 1;
	std::vector<Point> pnts;
	for(int i=p; i>=0; i--) {
//		std::cout << "Blue noise step " << i << "... " << std::flush;
		bool need_refinement;
		if(i == p) {
			// place initial points
			pnts = PlacePoints(mipmaps[i], i);
			need_refinement = true;
		}
		else {
			// split points
			pnts = Split(pnts, mipmaps[i], need_refinement);
		}
		// refine points for new density map
		if(need_refinement) {
			Refine(pnts, mipmaps[i], 2);
		}
//		std::cout << pnts.size() << " points." << std::endl;
		if(max_steps > 0 && p - i + 1 >= max_steps) {
			break;
		}
	}
	return pnts;
}

template<typename T>
void PlotPoints(const std::vector<Point>& points, const slimage::Image<T>& img, const slimage::Pixel<T>& color, bool plot_1px)
{
	for(Point p : points) {
		// round position
		int px = int(p.x + 0.5f);
		int py = int(p.y + 0.5f);
		if(px < 0 || int(img.width()) <= px || py < 0 || int(img.height()) <= py) {
			continue;
		}
		img(px, py) = color;
		if(!plot_1px) {
			// paint a star
			//    X
			//   XXX
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
			if(py + 1 < int(img.width())) {
				img(px, py+1) = color;
			}
		}
	}
}

void PlotPoints(const std::vector<Point>& points, const slimage::Image1ub& img, unsigned char grey, bool plot_1px)
{
	PlotPoints(points, img, slimage::Pixel1ub{grey}, plot_1px);
}

void PlotPoints(const std::vector<Point>& points, const slimage::Image3ub& img, const slimage::Pixel3ub& color, bool plot_1px)
{
	PlotPoints(points, img, color, plot_1px);
}

}

//----------------------------------------------------------------------------//

std::vector<Eigen::Vector2f> Fattal(const Eigen::MatrixXf& density)
{
	std::vector<fattal::Point> pnts = fattal::Compute(density, 0);
	std::vector<Eigen::Vector2f> v(pnts.size());
	std::transform(pnts.begin(), pnts.end(), v.begin(),
		[](const fattal::Point& p) {
			return Eigen::Vector2f(
				2.0f*static_cast<float>(p.x), 2.0f*static_cast<float>(p.y));
		});
	return v;
}

//----------------------------------------------------------------------------//
}
//----------------------------------------------------------------------------//
