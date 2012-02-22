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
