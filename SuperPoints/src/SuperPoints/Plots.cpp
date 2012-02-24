/*
 * Plotting.cpp
 *
 *  Created on: Feb 22, 2012
 *      Author: david
 */

#include "Superpixels.hpp"
#include "Plots.hpp"
#include <Danvil/SimpleEngine/Primitives.h>
#include <Danvil/SimpleEngine/GlHelpers.h>
#include <Slimage/Paint.hpp>
//----------------------------------------------------------------------------//
namespace dasp {
namespace plots {
//----------------------------------------------------------------------------//

void PlotClusterPoints(const slimage::Image3ub& img, const Cluster& cluster, const ImagePoints& points, const slimage::Pixel3ub& color)
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

void PlotClusters(const slimage::Image3ub& img, const Clustering& clustering, const std::vector<slimage::Pixel3ub>& colors)
{
	assert(clustering.cluster.size() == colors.size());
	for(unsigned int i=0; i<clustering.cluster.size(); i++) {
		PlotClusterPoints(img, clustering.cluster[i], clustering.points, colors[i]);
	}
}

void PlotClusterEllipse(const slimage::Image3ub& img, const Cluster& cluster, const slimage::Pixel3ub& color, bool filled)
{
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
	if(filled) {
		slimage::FillEllipse(img, cx, cy, p1x, p1y, p2x, p2y, color);
	}
	else {
		slimage::PaintEllipse(img, cx, cy, p1x, p1y, p2x, p2y, color);
	}

	//	Eigen::Vector2f p3 = opt.camera.project(cluster.center.world + 0.05f*cluster.center.normal);
	//	int p3x = std::round(p3[0]);
	//	int p3y = std::round(p3[1]);
	//	slimage::PaintLine(img, cx, cy, p3x, p3y, color);
}

void PlotEdges(const slimage::Image3ub& img, const std::vector<int>& labels, const slimage::Pixel3ub& color, unsigned int size)
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
			if(np > size) {
				img(x,y) = color;
				istaken[i] = true;
			}
		}
	}
}

template<typename K, unsigned int CC>
void PlotSeedsImpl(const slimage::Image<K,CC>& img, const std::vector<Seed>& seeds, const slimage::Pixel<K,CC>& color, int size)
{
	for(Seed s : seeds) {
		// round position
		slimage::PaintPoint(img, s.x, s.y, color, size);
	}
}

void PlotSeeds(const slimage::Image1ub& img, const std::vector<Seed>& seeds, unsigned char grey, int size)
{
	PlotSeedsImpl(img, seeds, slimage::Pixel1ub{grey}, size);
}

void PlotSeeds(const slimage::Image3ub& img, const std::vector<Seed>& seeds, const slimage::Pixel3ub& color, int size)
{
	PlotSeedsImpl(img, seeds, color, size);
}

namespace detail
{
	template<int M>
	slimage::Pixel3ub ComputePointColor(const Point&) {
		return {{0,0,0}};
	}

	template<>
	slimage::Pixel3ub ComputePointColor<UniBlack>(const Point& p) {
		return {{0,0,0}};
	}

	template<>
	slimage::Pixel3ub ComputePointColor<UniWhite>(const Point& p) {
		return {{255,255,255}};
	}

	template<>
	slimage::Pixel3ub ComputePointColor<Color>(const Point& p) {
		return RgbColor(p);
	}

	template<>
	slimage::Pixel3ub ComputePointColor<Depth>(const Point& p) {
		return DepthColor(p.depth_i16);
	}

	template<>
	slimage::Pixel3ub ComputePointColor<Gradient>(const Point& p) {
		if(p.isInvalid()) {
			return slimage::Pixel3ub{{0,0,0}};
		}
		return GradientColor(p.gradient);
	}

	template<int M>
	slimage::Pixel3ub ComputeClusterColor(const Cluster& c) {
		return ComputePointColor<M>(c.center);
	}

	template<>
	slimage::Pixel3ub ComputeClusterColor<Eccentricity>(const Cluster& c) {
		return IntensityColor(c.eccentricity, 0.0f, 1.0f);
	}

	template<>
	slimage::Pixel3ub ComputeClusterColor<Circularity>(const Cluster& c) {
		return IntensityColor(c.coverage, 0.0f, 1.0f);
	}

	template<>
	slimage::Pixel3ub ComputeClusterColor<Thickness>(const Cluster& c) {
		return IntensityColor(c.t, 0.0f, 0.01f);
	}

	template<int M>
	std::vector<slimage::Pixel3ub> ComputePixelColorsImpl(const dasp::ImagePoints& points)
	{
		std::vector<slimage::Pixel3ub> colors(points.size());
		for(unsigned int i=0; i<points.size(); i++) {
			colors[i] = ComputePointColor<M>(points[i]);
		}
		return colors;
	}

	template<int M>
	std::vector<slimage::Pixel3ub> ComputeClusterColorsImpl(const std::vector<Cluster>& clusters, const ClusterSelection& selection)
	{
		std::vector<slimage::Pixel3ub> colors(clusters.size());
		for(unsigned int i=0; i<clusters.size(); i++) {
			colors[i] = selection[i] ? ComputeClusterColor<M>(clusters[i]) : slimage::Pixel3ub{{0,0,0}};
		}
		return colors;
	}

	template<int M>
	void PlotPointsImpl(const slimage::Image3ub& img, const Clustering& c)
	{
		for(size_t i=0; i<img.getPixelCount(); i++) {
			img(i) = ComputePointColor<M>(c.points[i]);
		}
	}

}

#define ComputePixelColors_HELPER(T) case T: return detail::ComputePixelColorsImpl<T>(c.points);

std::vector<slimage::Pixel3ub> ComputePixelColors(const Clustering& c, ColorMode ccm)
{
	switch(ccm) {
	ComputePixelColors_HELPER(UniBlack)
	ComputePixelColors_HELPER(UniWhite)
	ComputePixelColors_HELPER(Color)
	ComputePixelColors_HELPER(Depth)
	ComputePixelColors_HELPER(Gradient)
	default: return detail::ComputePixelColorsImpl<-1>(c.points);
	}
}

#define ComputeClusterColors_HELPER(T) case T: return detail::ComputeClusterColorsImpl<T>(c.cluster, selection);

std::vector<slimage::Pixel3ub> ComputeClusterColors(const Clustering& c, ColorMode ccm, const ClusterSelection& selection)
{
	switch(ccm) {
	ComputeClusterColors_HELPER(UniBlack)
	ComputeClusterColors_HELPER(UniWhite)
	ComputeClusterColors_HELPER(Color)
	ComputeClusterColors_HELPER(Depth)
	ComputeClusterColors_HELPER(Gradient)
	ComputeClusterColors_HELPER(Eccentricity)
	ComputeClusterColors_HELPER(Circularity)
	ComputeClusterColors_HELPER(Thickness)
	default: return detail::ComputeClusterColorsImpl<-1>(c.cluster, selection);
	}
}

#define PlotPoints_HELPER(T) case T: detail::PlotPointsImpl<T>(img, c); break;

slimage::Image3ub PlotPoints(const Clustering& c, ColorMode cm)
{
	slimage::Image3ub img(c.width(), c.height());
	img.fill(0);
	switch(cm) {
	PlotPoints_HELPER(UniBlack)
	PlotPoints_HELPER(UniWhite)
	PlotPoints_HELPER(Color)
	PlotPoints_HELPER(Depth)
	PlotPoints_HELPER(Gradient)
	default: detail::PlotPointsImpl<-1>(img, c); break;
	}
	return img;
}

void PlotClusterCenters(const slimage::Image3ub& img, const Clustering& c, ColorMode ccm, int size, const ClusterSelection& selection)
{
	std::vector<slimage::Pixel3ub> colors = ComputeClusterColors(c, ccm, selection);
	for(size_t i=0; i<c.cluster.size(); i++) {
		if(selection[i]) {
			const Cluster& cluster = c.cluster[i];
			slimage::PaintPoint(img, cluster.center.spatial_x(), cluster.center.spatial_y(), colors[i], size);
		}
	}
}

void PlotClusterPoints(const slimage::Image3ub& img, const Clustering& c, ColorMode ccm, const ClusterSelection& selection)
{
	std::vector<slimage::Pixel3ub> colors = ComputeClusterColors(c, ccm, selection);
	for(size_t i=0; i<c.cluster.size(); i++) {
		if(selection[i]) {
			PlotClusterPoints(img, c.cluster[i], c.points, colors[i]);
		}
	}
}

void PlotClusterEllipses(const slimage::Image3ub& img, const Clustering& c, ColorMode ccm, const ClusterSelection& selection)
{
	std::vector<slimage::Pixel3ub> colors = ComputeClusterColors(c, ccm, selection);
	for(size_t i=0; i<c.cluster.size(); i++) {
		if(selection[i]) {
			PlotClusterEllipse(img, c.cluster[i], colors[i], false);
		}
	}
}

void PlotClusterEllipsesFilled(const slimage::Image3ub& img, const Clustering& c, ColorMode ccm, const ClusterSelection& selection)
{
	std::vector<slimage::Pixel3ub> colors = ComputeClusterColors(c, ccm, selection);
	for(size_t i=0; i<c.cluster.size(); i++) {
		if(selection[i]) {
			PlotClusterEllipse(img, c.cluster[i], colors[i], true);
		}
	}
}

void PlotClusters(slimage::Image3ub& img, const Clustering& c, ClusterMode mode, ColorMode cm, const ClusterSelection& selection)
{
	switch(mode) {
	case ClusterCenter: PlotClusterCenters(img, c, cm, 1, selection); break;
	default: case ClusterPoints: PlotClusterPoints(img, c, cm, selection); break;
	case ClusterEllipses: PlotClusterEllipses(img, c, cm, selection); break;
	case ClusterEllipsesFilled: PlotClusterEllipsesFilled(img, c, cm, selection); break;
	}
}

slimage::Image3ub PlotClusters(const Clustering& c, ClusterMode mode, ColorMode cm, const ClusterSelection& selection)
{
	slimage::Image3ub img(c.width(), c.height());
	img.fill(0);
	PlotClusters(img, c, mode, cm, selection);
	return img;
}

void RenderCluster(const Cluster& cluster, float r, const slimage::Pixel3ub& color)
{
	glDisable(GL_CULL_FACE);
	glPolygonMode(GL_FRONT, GL_FILL);
	glPolygonMode(GL_BACK, GL_LINE);
	// render circle at position
	Danvil::ctLinAlg::Vec3f pos = Danvil::ctLinAlg::Convert(cluster.center.world);
	Eigen::Vector3f v = cluster.center.world.cross(cluster.center.normal).normalized();
	Eigen::Vector3f u = v.cross(cluster.center.normal);
	Danvil::ctLinAlg::Vec3f major = r * Danvil::ctLinAlg::Convert(v);
	Danvil::ctLinAlg::Vec3f minor = r * Danvil::ctLinAlg::Convert(u);
	glColor3ub(color[0], color[1], color[2]);
	Danvil::SimpleEngine::Primitives::RenderEllipseCap(pos, major, minor);
	glColor3f(0.9f, 0.9f, 0.9f);
	Danvil::SimpleEngine::Primitives::RenderSegment(pos, pos + 0.02f * Danvil::ctLinAlg::Convert(cluster.center.normal));
}

void RenderClusters(const Clustering& clustering, ColorMode ccm, const ClusterSelection& selection)
{
	std::vector<slimage::Pixel3ub> colors = ComputeClusterColors(clustering, ccm, selection);
	for(unsigned int i=0; i<clustering.cluster.size(); i++) {
		if(selection[i]) {
			RenderCluster(clustering.cluster[i], clustering.opt.base_radius, colors[i]);
		}
	}
}

//----------------------------------------------------------------------------//
}}
//----------------------------------------------------------------------------//
