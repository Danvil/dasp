/*
 * BlueNoise.cpp
 *
 *  Created on: Feb 6, 2012
 *      Author: david
 */

#include "BlueNoise.hpp"
#include "Mipmaps.hpp"
#include <boost/random.hpp>
//----------------------------------------------------------------------------//
namespace dasp {
namespace BlueNoise {
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

float Energy(const std::vector<Point>& pnts, const slimage::Image1f& density)
{
	float error = 0.0f;
	for(unsigned int y=0; y<density.height(); y++) {
		for(unsigned int x=0; x<density.width(); x++) {
			float px = float(x);
			float py = float(y);
			float a = EnergyApproximation(pnts, px, py);
			float roh = density(x, y);
			error += std::abs(a - roh);
		}
	}
	return error;
}

void EnergyDerivative(const std::vector<Point>& pnts, const slimage::Image1f& density, unsigned int i, float& result_dE_x, float& result_dE_y)
{
	float dE_x = 0.0f;
	float dE_y = 0.0f;
	float px = pnts[i].x;
	float py = pnts[i].y;
	float ps = pnts[i].scale;
//		float ps_scl = 1.0f / ps;
	float ps_scl = 1.0f / (ps * ps);
	// find window (points outside the window do not affect the kernel)
	float radius = KernelRange * ps;
	float x_min = std::max(0, int(std::floor(px - radius)));
	float x_max = std::min(int(density.width()) - 1, int(std::ceil(px + radius)));
	float y_min = std::max(0, int(std::floor(py - radius)));
	float y_max = std::min(int(density.height()) - 1, int(std::ceil(py + radius)));
	// sum over window
	for(unsigned int y=y_min; y<=y_max; y++) {
		for(unsigned int x=x_min; x<=x_max; x++) {
			float ux = float(x);
			float uy = float(y);
			float dx = ux - px;
			float dy = uy - py;
//				float k_arg = std::sqrt(dx*dx + dy*dy) * ps_scl;
//				float k_val = KernelFunctor(k_arg);
			float k_arg_square = (dx*dx + dy*dy) * ps_scl;
			float k_val = KernelFunctorSquare(k_arg_square);
			float apx = EnergyApproximation(pnts, ux, uy);
			float roh = density(x, y);
			if(apx < roh) {
				k_val = -k_val;
			}
			dE_x += k_val * dx;
			dE_y += k_val * dy;
		}
	}
	float A = 1.0f / std::pow(ps, float(D + 1));
	result_dE_x = A * dE_x;
	result_dE_y = A * dE_y;
}

std::vector<Point> PlacePoints(const slimage::Image1f& density, unsigned int p)
{
	// access original index in a random order
	std::vector<unsigned int> indices(density.size());
	for(unsigned int i=0; i<indices.size(); i++) {
		indices[i] = i;
	}

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

	std::random_shuffle(indices.begin(), indices.end());
	// compute points
	std::vector<Point> pnts;
	pnts.reserve(indices.size());
	// compute current error in density
	float error_current = Energy(pnts, density);
	std::cout << "INITIAL ERROR: " << error_current << std::endl;
	// try add kernel points
	for(unsigned int i : indices) {
		float roh = density[i];
		if(roh == 0) {
//				std::cout << i << " roh is 0!" << std::endl;
			continue;
		}
		Point u;
		u.x = float(i % density.width());
		u.y = float(i / density.width());
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

void Refine(std::vector<Point>& points, const slimage::Image1f& density, unsigned int iterations)
{
	static boost::mt19937 rng;
	static boost::normal_distribution<float> rnd(0.0f, 1.0f); // standard normal distribution
	static boost::variate_generator<boost::mt19937&, boost::normal_distribution<float> > die(rng, rnd);
	constexpr float dt = 0.35f;
	constexpr float T = 0.0f;
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
			EnergyDerivative(points, density, i, dx, dy);
			p.x -= cA * dx;
			p.y -= cA * dy;
			if(T > 0.0f) {
				float cB = std::sqrt(T * c0);
				p.x += cB * die();
				p.y += cB * die();
			}
		}
	}
}

std::vector<Point> Split(const std::vector<Point>& points, const slimage::Image1f& density, bool& result_added)
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

std::vector<Point> Compute(const slimage::Image1f& density)
{
	// compute mipmaps
	std::vector<slimage::Image1f> mipmaps = Mipmaps::ComputeMipmaps(density, 16);
	int p = int(mipmaps.size()) - 1;
	std::vector<Point> pnts;
	for(int i=p; i>=0; i--) {
		std::cout << "Blue noise step " << i << "... " << std::flush;
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
			Refine(pnts, mipmaps[i], 16);
		}
		std::cout << pnts.size() << " points." << std::endl;
	}
	return pnts;
}

//----------------------------------------------------------------------------//
}}
//----------------------------------------------------------------------------//
