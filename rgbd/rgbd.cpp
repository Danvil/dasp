#include "rgbd.hpp"
#include <Slimage/IO.hpp>
#ifdef DASP_HAS_OPENNI
	#include "KinectGrabber.h"
#endif
#include <iostream>
#include <stdexcept>
#include <random>

constexpr int WIDTH = 640;
constexpr int HEIGHT = 480;

std::mt19937 random_engine;

slimage::Pixel3ub color(int r, int g, int b)
{
	return {{
		static_cast<unsigned char>(std::min(255, std::max(0, r))),
		static_cast<unsigned char>(std::min(255, std::max(0, g))),
		static_cast<unsigned char>(std::min(255, std::max(0, b)))}};
}

template<typename Dist>
slimage::Pixel3ub color_with_noise(int r, int g, int b, Dist d)
{
	return color(
		r + d(random_engine),
		g + d(random_engine),
		b + d(random_engine));
}

template<bool UseNoise=true>
class RgbdStreamTestUniform : public RgbdStream
{
public:
	RgbdStreamTestUniform() {
		data_.color = slimage::Image3ub(WIDTH, HEIGHT);
		data_.depth = slimage::Image1ui16(WIDTH, HEIGHT);
		constexpr int NC = UseNoise ? 15 : 0;
		constexpr int ND = UseNoise ? 50 : 0;
		std::uniform_int_distribution<int> dc(-NC,+NC);
		for(int i=0; i<HEIGHT; i++) {
			for(int j=0; j<WIDTH; j++) {
				int dd = static_cast<int>(
					ND*std::sin(30.0f*static_cast<float>(i)/static_cast<float>(HEIGHT))
					*std::sin(30.0f*static_cast<float>(j)/static_cast<float>(WIDTH))
				);
				data_.depth(j,i) = 1200 + dd;
				data_.color(j,i) = color_with_noise(32,128,128, dc);
			}
		}
	}
	bool grab() { return true; }
	Rgbd get() { return data_; }
private:
	Rgbd data_;
};

class RgbdStreamTestParaboloid : public RgbdStream
{
public:
	RgbdStreamTestParaboloid() {
		data_.color = slimage::Image3ub(WIDTH, HEIGHT, {{0,128,128}});
		data_.depth = slimage::Image1ui16(WIDTH, HEIGHT);
		for(int i=0; i<HEIGHT; i++) {
			int di = i - HEIGHT/2;
			for(int j=0; j<WIDTH; j++) {
				int dj = j - WIDTH/2;
				data_.depth(j,i) = static_cast<uint16_t>(600 + (di*di + dj*dj)/80);
			}
		}
	}
	bool grab() { return true; }
	Rgbd get() { return data_; }
private:
	Rgbd data_;
};

template<bool UseNoise=true>
class RgbdStreamTestSphere : public RgbdStream
{
public:
	RgbdStreamTestSphere() : time_(0) {}
	bool grab() { return true; }
	Rgbd get() {
		slimage::Image3ub color(WIDTH, HEIGHT);
		slimage::Image1ui16 depth(WIDTH, HEIGHT);
		constexpr int NC = UseNoise ? 15 : 0;
		constexpr int ND = UseNoise ? 50 : 0;
		std::uniform_int_distribution<int> dc(-NC,+NC);
		for(int i=0; i<HEIGHT; i++) {
			for(int j=0; j<WIDTH; j++) {
				float phi = 2.0f*3.1415f*static_cast<float>(time_)/100;
				const float r_move = 70.0f;
				int cx = WIDTH/2 + r_move*std::cos(phi);
				int cy = HEIGHT/2 + r_move*std::sin(phi);
				int dj = j - cx;
				int di = i - cy;
				int r = std::sqrt(di*di + dj*dj);
				int dd = static_cast<int>(
					ND*std::sin(30.0f*static_cast<float>(i)/static_cast<float>(HEIGHT))
					*std::sin(30.0f*static_cast<float>(j)/static_cast<float>(WIDTH))
				);
				if(r < 80) {
					depth(j,i) = 800 + dd;
					color(j,i) = color_with_noise(192,32,32, dc);
				}
				else {
					depth(j,i) = 1200 + dd;
					color(j,i) = color_with_noise(32,128,128, dc);
				}
			}
		}
		time_++;
		return {color, depth};
	}
private:
	int time_;
};

class RgbdStreamStatic : public RgbdStream
{
public:
	RgbdStreamStatic(const std::string& fn) {
		data_.color = slimage::Load3ub(fn + "_color.png");
		data_.depth = slimage::Load1ui16(fn + "_depth.pgm");
		if(data_.color.width() != data_.depth.width() || data_.color.height() != data_.depth.height()) {
			std::cerr << "WARNING: Size of color and depth image do not match!" << std::endl;
		}
	}
	bool grab() { return true; }
	Rgbd get() { return data_; }
private:
	Rgbd data_;
};

#ifdef DASP_HAS_OPENNI

class RgbdStreamOni : public RandomAccessRgbdStream
{
public:
	RgbdStreamOni(const std::string& fn) {
		grabber_ = std::make_shared<dasp::KinectGrabber>();
		grabber_->OpenFile(fn);
	}
	~RgbdStreamOni() {
		grabber_->Stop();
	}
	unsigned int numFrames() {
		return grabber_->NumFrames();
	}
	unsigned int tell() {
		return grabber_->TellFrame();
	}
	void seek(unsigned int frame) {
		grabber_->SeekToFrame(frame);
	}
	bool grab() {
		return grabber_->Grab();
	}
	Rgbd get() {
		return {
			grabber_->GetLastColor().clone(),
			grabber_->GetLastDepth().clone()
		};
	}
private:
	std::shared_ptr<dasp::KinectGrabber> grabber_;
};

class RgbdStreamKinectLive : public RgbdStream
{
public:
	RgbdStreamKinectLive(const std::string& fn_config) {
		grabber_ = std::make_shared<dasp::KinectGrabber>();
		grabber_->OpenConfig(fn_config);
	}
	~RgbdStreamKinectLive() {
		grabber_->Stop();
	}
	bool grab() {
		return grabber_->Grab();
	}
	Rgbd get() {
		return {
			grabber_->GetLastColor().clone(),
			grabber_->GetLastDepth().clone()
		};
	}
private:
	std::shared_ptr<dasp::KinectGrabber> grabber_;
};

#endif

std::shared_ptr<RgbdStream> FactorTest(const std::string& arg)
{
	if(arg == "uniform") {
		return std::make_shared<RgbdStreamTestUniform<true>>();
	}
	if(arg == "paraboloid") {
		return std::make_shared<RgbdStreamTestParaboloid>();
	}
	if(arg == "sphere") {
		return std::make_shared<RgbdStreamTestSphere<true>>();
	}
	std::cerr << "ERROR: Invalid rgbd test stream arg='" << arg << "'!" << std::endl;
	throw 0;
}

std::shared_ptr<RgbdStream> FactorStatic(const std::string& fn)
{
	return std::make_shared<RgbdStreamStatic>(fn);
}

std::shared_ptr<RandomAccessRgbdStream> FactorOni(const std::string& fn)
{
#ifdef DASP_HAS_OPENNI
	return std::make_shared<RgbdStreamOni>(fn);
#else
	std::cerr << "ERROR: Library configured without OpenNI! Enabel DASP_HAS_OPENNI in CMake." << std::endl;
	throw 0;
#endif
}

std::shared_ptr<RgbdStream> FactorKinectLive(const std::string& fn_config)
{
#ifdef DASP_HAS_OPENNI
	return std::make_shared<RgbdStreamKinectLive>(fn_config);
#else
	std::cerr << "ERROR: Library configured without OpenNI! Enabel DASP_HAS_OPENNI in CMake." << std::endl;
	throw 0;
#endif
}

std::shared_ptr<RgbdStream> FactorStream(const std::string& mode, const std::string& arg)
{
	if(mode == "test") {
		return FactorTest(arg);
	}
	if(mode == "static") {
		return FactorStatic(arg);
	}
	if(mode == "oni") {
		return FactorOni(arg);
	}
	if(mode == "live") {
		return FactorKinectLive(arg);
	}
	// default
	std::cerr << "ERROR: Invalid mode='" << mode << "' and arg='" << arg << "'!" << std::endl;
	throw 0;
}
