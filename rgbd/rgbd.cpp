#include "rgbd.hpp"
#include <Slimage/IO.hpp>
#ifdef DASP_HAS_OPENNI
	#include "KinectGrabber.h"
#endif
#include <iostream>

constexpr int WIDTH = 640;
constexpr int HEIGHT = 480;

class RgbdStreamTestUniform : public RgbdStream
{
public:
	RgbdStreamTestUniform() {
		data_.color = slimage::Image3ub(WIDTH, HEIGHT, {{0,128,128}});
		data_.depth = slimage::Image1ui16(WIDTH, HEIGHT);
		for(int i=0; i<HEIGHT; i++) {
			for(int j=0; j<WIDTH; j++) {
				data_.depth(j,i) = 1200;
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

class RgbdStreamTestSphere : public RgbdStream
{
public:
	RgbdStreamTestSphere() : time_(0) {}
	bool grab() { return true; }
	Rgbd get() {
		slimage::Image3ub color(WIDTH, HEIGHT, {{0,128,128}});
		slimage::Image1ui16 depth(WIDTH, HEIGHT);
		for(int i=0; i<HEIGHT; i++) {
			for(int j=0; j<WIDTH; j++) {
				depth(j,i) = 1200;
				float phi = 2.0f*3.1415f*static_cast<float>(time_)/100;
				const float r_move = 70.0f;
				int cx = WIDTH/2 + r_move*std::cos(phi);
				int cy = HEIGHT/2 + r_move*std::sin(phi);
				int dj = j - cx;
				int di = i - cy;
				int d = std::sqrt(di*di + dj*dj);
				if(d < 80) {
					depth(j,i) = 800;
					color(j,i) = {{255,0,0}};
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

std::shared_ptr<RgbdStream> FactorTest(const std::string& tag)
{
	if(tag == "uniform") {
		return std::make_shared<RgbdStreamTestUniform>();
	}
	if(tag == "paraboloid") {
		return std::make_shared<RgbdStreamTestParaboloid>();
	}
	if(tag == "sphere") {
		return std::make_shared<RgbdStreamTestSphere>();
	}
	return std::shared_ptr<RgbdStream>();
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
	return std::shared_ptr<RandomAccessRgbdStream>();
#endif
}

std::shared_ptr<RgbdStream> FactorKinectLive(const std::string& fn_config)
{
#ifdef DASP_HAS_OPENNI
	return std::make_shared<RgbdStreamKinectLive>(fn_config);
#else
	return std::shared_ptr<RgbdStream>();
#endif
}

