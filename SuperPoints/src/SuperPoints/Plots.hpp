/*
 * Plotting.hpp
 *
 *  Created on: Feb 22, 2012
 *      Author: david
 */

#ifndef PLOTTING_HPP_
#define PLOTTING_HPP_
//----------------------------------------------------------------------------//
#include "Superpixels.hpp"
#include <Danvil/Color.h>
//----------------------------------------------------------------------------//
namespace dasp {
namespace plots {
//----------------------------------------------------------------------------//

/** Get point color */
inline slimage::Pixel3ub RgbColor(const Point& p) {
	return {{
		static_cast<unsigned char>(std::min(255.0f, std::max(0.0f, 255.0f*p.color[0]))),
		static_cast<unsigned char>(std::min(255.0f, std::max(0.0f, 255.0f*p.color[1]))),
		static_cast<unsigned char>(std::min(255.0f, std::max(0.0f, 255.0f*p.color[2])))
	}};
}

/** Color visualization of gradient */
inline slimage::Pixel3ub GradientColor(const Eigen::Vector2f& g)
{
	float x = std::max(0.0f, std::min(1.0f, 0.5f + g[0]));
	float y = std::max(0.0f, std::min(1.0f, 0.5f + g[1]));
	return slimage::Pixel3ub{{
			static_cast<unsigned char>(255.0f*0.5f*(1.0f - x + y)),
			static_cast<unsigned char>(255.0f*0.5f*(2.0f - x - y)),
			static_cast<unsigned char>(255.0f*0.5f*(x + y))}};
}

/** Color visualization of kinect depth */
inline slimage::Pixel3ub DepthColor(uint16_t d16)
{
	// base gradient: blue -> red -> yellow
	static auto cm = Danvil::ContinuousIntervalColorMapping<unsigned char, uint16_t>::Factor_Blue_Red_Yellow();
	cm.setRange(400,2000);
	if(d16 == 0) {
		return slimage::Pixel3ub{{0,0,0}};
	}
	else {
		Danvil::Colorub color = cm(d16);
		unsigned int q = d16 % 25;
		unsigned char r = std::max(0, int(color.r) - int(q));
		unsigned char g = std::max(0, int(color.g) - int(q));
		unsigned char b = std::max(0, int(color.b) - int(q));
		return slimage::Pixel3ub{{r,g,b}};
	}
}

/** Color visualization of intensity */
inline slimage::Pixel3ub IntensityColor(float x, float min=0.0f, float max=1.0f)
{
	static auto cm = Danvil::ContinuousIntervalColorMapping<unsigned char, float>::Factor_Blue_Red_Yellow();
	cm.setRange(min, max);
	Danvil::Colorub color = cm(x);
	return slimage::Pixel3ub{{color.r,color.g,color.b}};
}

void PlotClusterPoints(const slimage::Image3ub& img, const Cluster& cluster, const ImagePoints& points, const slimage::Pixel3ub& color);

void PlotClusters(const slimage::Image3ub& img, const Clustering& clustering, const std::vector<slimage::Pixel3ub>& colors);

void PlotClusterEllipse(const slimage::Image3ub& img, const Cluster& cluster, const slimage::Pixel3ub& color, bool filled);

void PlotEdges(const slimage::Image3ub& img, const std::vector<int>& point_labels, const slimage::Pixel3ub& color, unsigned int size=1);

void PlotSeeds(const slimage::Image1ub& img, const std::vector<Seed>& seeds, unsigned char grey=0, int size=1);

void PlotSeeds(const slimage::Image3ub& img, const std::vector<Seed>& seeds, const slimage::Pixel3ub& color=slimage::Pixel3ub{{0,0,0}}, int size=1);

enum ColorMode {
	UniBlack,
	UniWhite,
	Color,
	Depth,
	Gradient,
	Eccentricity,
	Circularity,
	Thickness
};

std::vector<slimage::Pixel3ub> ComputePixelColors(const Clustering& c, ColorMode ccm);

std::vector<slimage::Pixel3ub> ComputeClusterColors(const Clustering& c, ColorMode ccm);

slimage::Image3ub PlotPoints(const Clustering& c, ColorMode pcm);

void PlotClusterCenters(const slimage::Image3ub& img, const Clustering& c, ColorMode ccm, int size);

void PlotClusterPoints(const slimage::Image3ub& img, const Clustering& c, ColorMode ccm);

void PlotClusterEllipses(const slimage::Image3ub& img, const Clustering& c, ColorMode ccm);

void PlotClusterEllipsesFilled(const slimage::Image3ub& img, const Clustering& c, ColorMode ccm);

enum ClusterMode {
	ClusterCenter,
	ClusterPoints,
	ClusterEllipses,
	ClusterEllipsesFilled
};

void PlotClusters(slimage::Image3ub& img, const Clustering& c, ClusterMode cm, ColorMode ccm);

slimage::Image3ub PlotClusters(const Clustering& c, ClusterMode cm, ColorMode ccm);

//----------------------------------------------------------------------------//
}}
//----------------------------------------------------------------------------//
#endif
