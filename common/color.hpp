#ifndef COMMON_COLOR_HPP_
#define COMMON_COLOR_HPP_
//----------------------------------------------------------------------------//
#include <Danvil/Color.h>
#include <Eigen/Dense>
#include <Slimage/Slimage.hpp>
//----------------------------------------------------------------------------//
namespace common {
//----------------------------------------------------------------------------//

template<typename K>
inline Eigen::Vector3f GreyColor(K x, K a, K b)
{
	float p = (static_cast<float>(x) - static_cast<float>(a)) / (static_cast<float>(b) - static_cast<float>(a));
	p = std::min(1.0f, std::max(0.0f, p));
	return {p,p,p};
}

/** Computes a color to express a similarity value between 0 and 1 */
inline Eigen::Vector3f SimilarityColor(float x)
{
	static auto cm = Danvil::ContinuousIntervalColorMapping<float, float>::Factor_Black_Blue_Red_Yellow_White();
	cm.setRange(0.0f, 1.0f);
	Danvil::Colorf color = cm(x);
	return {color.r,color.g,color.b};
}

inline Eigen::Vector3f IntensityColor(float x, float a=0.0f, float b=1.0f)
{
	static auto cm = Danvil::ContinuousIntervalColorMapping<float, float>::Factor_Black_Blue_Red_Yellow_White();
	cm.setRange(a, b);
	Danvil::Colorf color = cm(x);
	return {color.r,color.g,color.b};
}

inline Eigen::Vector3f PlusMinusColor(float x, float range=1.0f)
{
	static auto cm = Danvil::ContinuousIntervalColorMapping<float, float>::Factor_MinusPlus();
	cm.setRange(-range, +range);
	Danvil::Colorf color = cm(x);
	return {color.r,color.g,color.b};
}

template<typename T>
inline Eigen::Vector3f CountColor(T num, T min, T max)
{
	return IntensityColor(
		(num < min)
			? 0.0f
			: static_cast<float>(num - min)/static_cast<float>(max-min)
	);
}

inline slimage::Pixel3ub ColorToPixel(const Eigen::Vector3f& color) {
	return slimage::Pixel3ub{{
		static_cast<unsigned char>(255.f*std::min(1.0f, std::max(0.0f, color[0]))),
		static_cast<unsigned char>(255.f*std::min(1.0f, std::max(0.0f, color[1]))),
		static_cast<unsigned char>(255.f*std::min(1.0f, std::max(0.0f, color[2])))
	}};
}

inline slimage::Image3ub ColorizeDepth(const slimage::Image1ui16& img16, uint16_t min, uint16_t max) {
	slimage::Image3ub img(img16.width(), img16.height());
	const int n = img16.size();
	for(int i=0; i<n; i++) {
		img[i] = ColorToPixel(CountColor((uint16_t)img16[i], min, max));
	}
	return img;
}

inline slimage::Image3ub GreyDepth(const slimage::Image1ui16& img16, uint16_t min, uint16_t max) {
	slimage::Image3ub img(img16.width(), img16.height());
	const int n = img16.size();
	for(int i=0; i<n; i++) {
		img[i] = ColorToPixel(GreyColor<uint16_t>(img16[i], min, max));
	}
	return img;
}

template<typename CF>
inline slimage::Image3ub MatrixToImage(const Eigen::MatrixXf& mat, CF cf)
{
	slimage::Image3ub vis = slimage::Image3ub(mat.rows(), mat.cols());
	const float* p = mat.data();
	for(unsigned int i=0; i<vis.size(); i++) {
		vis[i] = ColorToPixel(cf(p[i]));
	}
	return vis;
}

//----------------------------------------------------------------------------//
}
//----------------------------------------------------------------------------//
#endif
