/*
 * BlueNoise.hpp
 *
 *  Created on: Feb 6, 2012
 *      Author: david
 */

#ifndef BLUENOISE_HPP_
#define BLUENOISE_HPP_
//----------------------------------------------------------------------------//
#include <Slimage/Slimage.hpp>
#include <Danvil/Tools/FunctionCache.h>
#include <vector>
#include <cmath>
//----------------------------------------------------------------------------//
namespace dasp {
//----------------------------------------------------------------------------//

namespace Quadtree
{

	template<typename T>
	struct Element
	{
		float x, y;
		T data;
	};

	/**
	 * 0 1
	 * 2 3
	 */
	template<typename T>
	struct Node
	{
		Node(float x, float y, float s)
		: center_x_(x), center_y_(y), size_(s), is_leaf_(true) {}

		template<typename F>
		void for_each_in_range(float x, float y, F f) {
			// FIXME
		}

		template<typename F>
		void for_each(F f) {
			if(is_leaf_) {
				std::for_each(data_.begin(), data_.end(), [f](const Element<T>& e) { f(e.data); } );
			}
			else {
				for(unsigned int i=0; i<4; i++) {
					children_[i]->for_each(f);
				}
			}
		}

		void expand(unsigned int lvl=1) {
			if(is_leaf_) {
				float a = 0.25f * size_;
				float s = 0.50f * size_;
				children_[0].reset(new Node(center_x_ - a, center_y_ - a, s));
				children_[1].reset(new Node(center_x_ + a, center_y_ - a, s));
				children_[2].reset(new Node(center_x_ - a, center_y_ + a, s));
				children_[3].reset(new Node(center_x_ + a, center_y_ + a, s));
				// assign data
				for(const Element<T>& e : data_) {
					children_[findQuadrant(e)]->add(e);
				}
				data_.clear();
				is_leaf_ = false;
			}
			if(lvl > 1) {
				for(unsigned int i=0; i<4; i++) {
					children_[i]->expand(lvl - 1);
				}
			}
		}

		void findQuadrant(const Element<T>& e) {
			bool tx = (e.x < center_x_);
			bool ty = (e.y < center_y_);
			return ty
					? (tx ? 0 : 1)
					: (tx ? 2 : 3);
		}

		void add(const std::vector<Element<T>>& ve) {
			for(const std::vector<Element<T>>& e : ve) {
				add(e);
			}
		}

		void add(const Element<T>& e) {
			if(is_leaf_) {
				data_.push_back(e);
			}
			else {
				children_[findQuadrant(e)]->add(e);
			}
		}

	private:
		float center_x_, center_y_;
		float size_;
		bool is_leaf_;
		boost::shared_ptr<Node> children_[4];
		std::vector<Element<T>> data_;
	};

}

namespace BlueNoise
{
	// need to change some other functions too!!!
	constexpr unsigned int D = 2;

	struct Point {
		float x, y;
		float weight;
		float scale;
	};

	constexpr float KernelRange = 2.5f;

	constexpr float cMaxRefinementScale = 10.0f;

	constexpr float cPi = 3.141592654f;

	/** phi(x) = exp(-pi*x*x) */
	inline
	float KernelFunctorImpl(float d) {
		return std::exp(-cPi*d*d);
	}

	inline
	float KernelFunctor(float d) {
		static Danvil::FunctionCache<float,1> cache(0.0f, KernelRange, &KernelFunctorImpl);
		return cache(std::abs(d));
	}

	inline
	float KernelFunctorSquareImpl(float d) {
		return std::exp(-cPi*d);
	}

	inline
	float KernelFunctorSquare(float d) {
		static Danvil::FunctionCache<float,1> cache(0.0f, KernelRange*KernelRange, &KernelFunctorSquareImpl);
		return cache(d);
	}

	inline
	float ZeroBorderAccess(const slimage::Image1f& density, int x, int y) {
		if(0 <= x && x < int(density.width()) && 0 <= y && y < int(density.height())) {
			return density(x, y);
		}
		else {
			return 0.0f;
		}
	}

	inline
	float KernelScaleFunction(float roh, float weight) {
//		return std::pow(roh / weight, -1.0f / float(D));
		return 1.0f / std::sqrt(roh / weight);
	}

	inline
	float ScalePowerD(float s) {
		//return std::pow(s, -float(D));
		return 1.0f / (s*s);
	}

	float EnergyApproximation(const Quadtree::Node<Point>& pnts, float x, float y);

	float Energy(const Quadtree::Node<Point>& pnts, const slimage::Image1f& density);

	float EnergyDerivative(const Quadtree::Node<Point>& pnts, const slimage::Image1f& density, unsigned int i, float& result_dE_x, float& result_dE_y);

	std::vector<Point> PlacePoints(const slimage::Image1f& density, unsigned int p);

	void Refine(std::vector<Point>& points, const slimage::Image1f& density, unsigned int iterations);

	std::vector<Point> Split(const std::vector<Point>& points, const slimage::Image1f& density, bool& result_added);

	std::vector<Point> Compute(const slimage::Image1f& density, unsigned int max_steps=0);

	void PlotPoints(const std::vector<Point>& points, const slimage::Image1ub& img, bool plot_1px=true);

}

//----------------------------------------------------------------------------//
}
//----------------------------------------------------------------------------//
#endif
