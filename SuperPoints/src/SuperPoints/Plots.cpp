/*
 * Plotting.cpp
 *
 *  Created on: Feb 22, 2012
 *      Author: david
 */

#include "Superpixels.hpp"
#include "Plots.hpp"
#include <Slimage/Paint.hpp>
#include <Danvil/Color.h>
//----------------------------------------------------------------------------//
namespace dasp {
namespace plots {
//----------------------------------------------------------------------------//

namespace detail
{
	slimage::Pixel3ub GradientColor(const Eigen::Vector2f& g)
	{
		float x = std::max(0.0f, std::min(1.0f, 0.5f + g[0]));
		float y = std::max(0.0f, std::min(1.0f, 0.5f + g[1]));
		return slimage::Pixel3ub{{
				static_cast<unsigned char>(255.0f*0.5f*(1.0f - x + y)),
				static_cast<unsigned char>(255.0f*0.5f*(2.0f - x - y)),
				static_cast<unsigned char>(255.0f*0.5f*(x + y))}};
	}

	slimage::Pixel3ub DepthColor(uint16_t d16)
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

	slimage::Pixel3ub IntensityColor(float x, float min=0.0f, float max=1.0f)
	{
		// base gradient: blue -> red -> yellow
		static auto cm = Danvil::ContinuousIntervalColorMapping<unsigned char, float>::Factor_Blue_Red_Yellow();
		cm.setRange(min, max);
		Danvil::Colorub color = cm(x);
		return slimage::Pixel3ub{{color.r,color.g,color.b}};
	}
}

void PlotCluster(const Cluster& cluster, const ImagePoints& points, const slimage::Image3ub& img)
{
	unsigned char c_col_r = 255.0f * cluster.center.color[0];
	unsigned char c_col_g = 255.0f * cluster.center.color[1];
	unsigned char c_col_b = 255.0f * cluster.center.color[2];
	PlotCluster(cluster, points, img, slimage::Pixel3ub{{c_col_r, c_col_g, c_col_b}});
}

void PlotCluster(const Cluster& cluster, const ImagePoints& points, const slimage::Image3ub& img, const slimage::Pixel3ub& color)
{
	assert(cluster.hasPoints());
	// plot all pixels belonging to the cluster in the color of the cluster center
	for(unsigned int i : cluster.pixel_ids) {
		const Point& p = points[i];
		img(p.spatial_x(), p.spatial_y()) = color;
	}
	// plot the cluster center (using some kind of inverse color)
	int cx = cluster.center.spatial_x();
	int cy = cluster.center.spatial_y();
	if(0 <= cx && cx < int(points.width()) && 0 <= cy && cy < int(points.height())) {
		unsigned char* col = img.pointer(cx, cy);
		col[0] = 255 - color[0];
		col[1] = 255 - color[1];
		col[2] = 255 - color[2];
	}

}

void PlotCluster(const std::vector<Cluster>& clusters, const ImagePoints& points, const slimage::Image3ub& img)
{
//	std::vector<Danvil::ColorUB> cols = {
//			Danvil::Color::Red, Danvil::Color::Green, Danvil::Color::Blue, Danvil::Color::Yellow, Danvil::Color::Cyan, Danvil::Color::Magenta, Danvil::Color::Black
//	};
	img.fill(0);
	for(const Cluster& x : clusters) {
		PlotCluster(x, points, img);
	}
}

void PlotClusterCross(const Cluster& cluster, const slimage::Image3ub& img, const Parameters& opt)
{
	unsigned char c_col_r = 255.0f * cluster.center.color[0];
	unsigned char c_col_g = 255.0f * cluster.center.color[1];
	unsigned char c_col_b = 255.0f * cluster.center.color[2];
	slimage::Pixel3ub color{{c_col_r,c_col_g,c_col_b}};

	int cx = cluster.center.spatial_x();
	int cy = cluster.center.spatial_y();

	Eigen::Vector2f g0 = cluster.center.gradient;
	float g_norm_sq = g0.squaredNorm();
	if(g_norm_sq == 0.0f) {
		g0 = Eigen::Vector2f::Unit(0);
	}
	else {
		g0 /= std::sqrt(g_norm_sq);
	}
	float sp_0 = cluster.center.image_super_radius;
	float sp_small = sp_0 / std::sqrt(g_norm_sq + 1.0f);
	int p1x = static_cast<int>(sp_small * g0[0]);
	int p1y = static_cast<int>(sp_small * g0[1]);
	int p2x = static_cast<int>(- sp_0 * g0[1]);
	int p2y = static_cast<int>(+ sp_0 * g0[0]);
//	std::cout << cluster.center.gradient.transpose() << " --- " << sp_0 << " --- " << sp_small << std::endl;
//	std::cout << p1x << "," << p1y << " --- " << p2x << "," << p2y << std::endl;
//	slimage::PaintLine(img, cx, cy, cx + p1x, cy + p1y, color);
//	slimage::PaintLine(img, cx, cy, cx + p2x, cy + p2y, color);
	slimage::PaintEllipse(img, cx, cy, p1x, p1y, p2x, p2y, color);

	//	Eigen::Vector2f p3 = opt.camera.project(cluster.center.world + 0.05f*cluster.center.normal);
	//	int p3x = std::round(p3[0]);
	//	int p3y = std::round(p3[1]);
	//	slimage::PaintLine(img, cx, cy, p3x, p3y, color);
}

void PlotClustersCross(const std::vector<Cluster>& clusters, const slimage::Image3ub& img, const Parameters& opt)
{
	for(const Cluster& x : clusters) {
		PlotClusterCross(x, img, opt);
	}
}

void PlotEdges(const std::vector<int>& labels, const slimage::Image3ub& img, unsigned int edge_w, unsigned char edge_r, unsigned char edge_g, unsigned char edge_b)
{
	const int dx8[8] = {-1, -1,  0,  1, 1, 1, 0, -1};
	const int dy8[8] = { 0, -1, -1, -1, 0, 1, 1,  1};

	unsigned int width = img.width();
	unsigned int height = img.height();

	std::vector<bool> istaken(width*height, false);

	unsigned int i = 0;
	for(unsigned int y=0; y<height; y++) {
		for(unsigned int x=0; x<width; x++, i++) {
			unsigned int np = 0;
			for(unsigned int k=0; k<8; k++) {
				int sx = int(x) + dx8[k];
				int sy = int(y) + dy8[k];
				if(sx >= 0 && sx < int(width) && sy >= 0 && sy < int(height)) {
					int si = sy*width + sx;
					//if( false == istaken[si] )//comment this to obtain internal contours
					{
						if(labels[i] != labels[si]) {
							np++;
						}
					}
				}
			}
			if(np > edge_w) {
				img(x,y) = slimage::Pixel3ub{{edge_r, edge_g, edge_b}};
				istaken[i] = true;
			}
		}
	}
}

template<typename K, unsigned int CC>
void PlotSeedsImpl(const std::vector<Seed>& seeds, const slimage::Image<K,CC>& img, const slimage::Pixel<K,CC>& color, int size)
{
	for(Seed s : seeds) {
		// round position
		int px = s.x;
		int py = s.y;
		if(px < 0 || int(img.width()) <= px || py < 0 || int(img.height()) <= py) {
			continue;
		}
		if(size == 1) {
			img(px, py) = color;

		}
		else if(size == 2) {
			// paint a star
			//    X
			//   X X
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
			if(py + 1 < int(img.height())) {
				img(px, py+1) = color;
			}
		}
		else {
			// paint a circle
			//    X
			//   X X
			//  X   X
			//   X X
			//    X
			if(1 <= px && 1 <= py) {
				img(px-1, py-1) = color;
			}
			if(1 <= px && py + 1 < int(img.height())) {
				img(px-1, py+1) = color;
			}
			if(px + 1 < int(img.width()) && 1 <= py) {
				img(px+1, py-1) = color;
			}
			if(px + 1 < int(img.width()) && py + 1 < int(img.height())) {
				img(px+1, py+1) = color;
			}
			if(2 <= px) {
				img(px-2, py) = color;
			}
			if(px + 2 < int(img.width())) {
				img(px+2, py) = color;
			}
			if(2 <= py) {
				img(px, py-2) = color;
			}
			if(py + 2 < int(img.height())) {
				img(px, py+2) = color;
			}
		}
	}
}

void PlotSeeds(const std::vector<Seed>& seeds, const slimage::Image1ub& img, unsigned char grey, int size)
{
	PlotSeedsImpl(seeds, img, slimage::Pixel1ub{grey}, size);
}

void PlotSeeds(const std::vector<Seed>& seeds, const slimage::Image3ub& img, const slimage::Pixel3ub& color, int size)
{
	PlotSeedsImpl(seeds, img, color, size);
}

slimage::Image3ub PlotPoints(const Clustering& c, ColorMode pcm)
{
	std::vector<slimage::Pixel3ub> colors = ComputePixelColors(c, pcm);
	slimage::Image3ub img(c.points.width(), c.points.height());
	// FIXME implement
	return img;
}

slimage::Image3ub PlotClusters(const Clustering& c, ClusterMode cm, ColorMode ccm)
{
	std::vector<slimage::Pixel3ub> colors = ComputeClusterColors(c, ccm);
	slimage::Image3ub img(c.points.width(), c.points.height());
	// FIXME implement
	return img;
}

namespace detail
{
	template<int M>
	slimage::Pixel3ub ComputePointColor(const Point&) {
		return {{0,0,0}};
	}

	template<>
	slimage::Pixel3ub ComputePointColor<Color>(const Point& p) {
		return {{
			static_cast<unsigned char>(std::min(255.0f, std::max(0.0f, 255.0f*p.color[0]))),
			static_cast<unsigned char>(std::min(255.0f, std::max(0.0f, 255.0f*p.color[1]))),
			static_cast<unsigned char>(std::min(255.0f, std::max(0.0f, 255.0f*p.color[2])))
		}};
	}

	template<>
	slimage::Pixel3ub ComputePointColor<Depth>(const Point& p) {
		return DepthColor(p.depth_i16);
	}

	template<>
	slimage::Pixel3ub ComputePointColor<Gradient>(const Point& p) {
		return GradientColor(p.gradient);
	}

	template<int M>
	std::vector<slimage::Pixel3ub> ComputePixelColors(const dasp::ImagePoints& points, ColorMode ccm)
	{
		std::vector<slimage::Pixel3ub> colors(points.size());
		for(unsigned int i=0; i<points.size(); i++) {
			colors[i] = ComputePointColor<M>(points[i]);
		}
		return colors;
	}

	template<int M>
	slimage::Pixel3ub ComputeClusterColor(const Cluster& c) {
		return ComputePointColor<M>(c.center);
	}

//	template<>
//	slimage::Pixel3ub ComputeClusterColor<Eccentricity>(const Cluster& c) {
//		return IntensityColor(c.)
//	}

	template<int M>
	std::vector<slimage::Pixel3ub> ComputeClusterColors(const std::vector<Cluster>& clusters, ColorMode ccm)
	{
		std::vector<slimage::Pixel3ub> colors(clusters.size());
		for(unsigned int i=0; i<clusters.size(); i++) {
			colors[i] = ComputeClusterColor<M>(clusters[i]);
		}
		return colors;
	}

}

#define ComputePixelColors_HELPER(T) case T: return detail::ComputePixelColors<T>(c.points, ccm);

std::vector<slimage::Pixel3ub> ComputePixelColors(const Clustering& c, ColorMode ccm)
{
	switch(ccm) {
	ComputePixelColors_HELPER(Color)
	ComputePixelColors_HELPER(Depth)
	ComputePixelColors_HELPER(Gradient)
	default: return detail::ComputePixelColors<-1>(c.points, ccm);
	}
}

#define ComputeClusterColors_HELPER(T) case T: return detail::ComputeClusterColors<T>(c.cluster, ccm);

std::vector<slimage::Pixel3ub> ComputeClusterColors(const Clustering& c, ColorMode ccm)
{
	switch(ccm) {
	ComputeClusterColors_HELPER(Color)
	ComputeClusterColors_HELPER(Depth)
	ComputeClusterColors_HELPER(Gradient)
	ComputeClusterColors_HELPER(Eccentricity)
	ComputeClusterColors_HELPER(Circularity)
	ComputeClusterColors_HELPER(Thickness)
	default: return detail::ComputeClusterColors<-1>(c.cluster, ccm);
	}
}

//----------------------------------------------------------------------------//
}}
//----------------------------------------------------------------------------//
