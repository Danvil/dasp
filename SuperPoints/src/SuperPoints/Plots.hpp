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

inline slimage::Pixel3ub RgbColor(const Point& p) {
	return {{
		static_cast<unsigned char>(std::min(255.0f, std::max(0.0f, 255.0f*p.color[0]))),
		static_cast<unsigned char>(std::min(255.0f, std::max(0.0f, 255.0f*p.color[1]))),
		static_cast<unsigned char>(std::min(255.0f, std::max(0.0f, 255.0f*p.color[2])))
	}};
}

inline slimage::Pixel3ub GradientColor(const Eigen::Vector2f& g)
{
	float x = std::max(0.0f, std::min(1.0f, 0.5f + g[0]));
	float y = std::max(0.0f, std::min(1.0f, 0.5f + g[1]));
	return slimage::Pixel3ub{{
			static_cast<unsigned char>(255.0f*0.5f*(1.0f - x + y)),
			static_cast<unsigned char>(255.0f*0.5f*(2.0f - x - y)),
			static_cast<unsigned char>(255.0f*0.5f*(x + y))}};
}

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

inline slimage::Pixel3ub IntensityColor(float x, float min=0.0f, float max=1.0f)
{
	// base gradient: blue -> red -> yellow
	static auto cm = Danvil::ContinuousIntervalColorMapping<unsigned char, float>::Factor_Blue_Red_Yellow();
	cm.setRange(min, max);
	Danvil::Colorub color = cm(x);
	return slimage::Pixel3ub{{color.r,color.g,color.b}};
}

void PlotCluster(const Cluster& cluster, const ImagePoints& points, const slimage::Image3ub& img, const slimage::Pixel3ub& color);

void PlotCluster(const Clustering& clustering, const slimage::Image3ub& img, const std::vector<slimage::Pixel3ub>& colors);

void PlotClusterCross(const Cluster& cluster, const slimage::Image3ub& img, const Parameters& opt);

void PlotClustersCross(const Clustering& clustering, const slimage::Image3ub& img);

void PlotEdges(const std::vector<int>& point_labels, const slimage::Image3ub& img, unsigned int edge_w, unsigned char edge_r, unsigned char edge_g, unsigned char edge_b);

void PlotSeeds(const std::vector<Seed>& seeds, const slimage::Image1ub& img, unsigned char grey=0, int size=1);

void PlotSeeds(const std::vector<Seed>& seeds, const slimage::Image3ub& img, const slimage::Pixel3ub& color=slimage::Pixel3ub{{0,0,0}}, int size=1);

enum ColorMode {
	Color,
	Depth,
	Gradient,
	Eccentricity,
	Circularity,
	Thickness
};

enum ClusterMode {
	ClusterCenter,
	ClusterPixels,
	ClusterEllipses,
	ClusterEllipsesFilled
};

slimage::Image3ub PlotPoints(const Clustering& c, ColorMode pcm);

slimage::Image3ub PlotClusters(const Clustering& c, ClusterMode cm, ColorMode ccm);

std::vector<slimage::Pixel3ub> ComputePixelColors(const Clustering& c, ColorMode ccm);

std::vector<slimage::Pixel3ub> ComputeClusterColors(const Clustering& c, ColorMode ccm);

//----------------------------------------------------------------------------//
}}
//----------------------------------------------------------------------------//
#endif
