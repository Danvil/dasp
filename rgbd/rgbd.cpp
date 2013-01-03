#include "rgbd.hpp"
#include <Slimage/IO.hpp>
#include <kinect/KinectGrabber.h>

constexpr int WIDTH = 640;
constexpr int HEIGHT = 480;

class ScenarioTestUniform : public Scenario
{
public:
	ScenarioTestUniform() {
		data_.color = slimage::Image3ub(WIDTH, HEIGHT, {{0,128,128}});
		data_.depth = slimage::Image1ui16(WIDTH, HEIGHT);
		for(int i=0; i<HEIGHT; i++) {
			for(int j=0; j<WIDTH; j++) {
				data_.depth(j,i) = 1200;
			}
		}
	}
	Rgbd pop() { return data_; }
private:
	Rgbd data_;
};

class ScenarioTestParaboloid : public Scenario
{
public:
	ScenarioTestParaboloid() {
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
	Rgbd pop() { return data_; }
private:
	Rgbd data_;
};

class ScenarioTestSphere : public Scenario
{
public:
	ScenarioTestSphere() : time_(0) {}
	Rgbd pop() {
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

class ScenarioStatic : public Scenario
{
public:
	ScenarioStatic(const std::string& fn) {
		data_.color = slimage::Load3ub(fn + "_color.png");
		data_.depth = slimage::Load1ui16(fn + "_depth.pgm");
	}
	Rgbd pop() { return data_; }
private:
	Rgbd data_;
};

class ScenarioOni : public Scenario
{
public:
	ScenarioOni(const std::string& fn, int offset=0) {
		grabber_ = std::make_shared<dasp::KinectGrabber>();
		grabber_->OpenFile(fn);
		grabber_->SeekToFrame(offset);
	}
	Rgbd pop() {
		Rgbd data;
		bool ok = grabber_->Grab();
		if(ok) {
			data.color = grabber_->GetLastColor().clone();
			data.depth = grabber_->GetLastDepth().clone();
		}
		return data;
	}
private:
	std::shared_ptr<dasp::KinectGrabber> grabber_;
};

std::shared_ptr<Scenario> FactorTest(const std::string& tag)
{
	if(tag == "uniform") {
		return std::make_shared<ScenarioTestUniform>();
	}
	if(tag == "paraboloid") {
		return std::make_shared<ScenarioTestParaboloid>();
	}
	if(tag == "sphere") {
		return std::make_shared<ScenarioTestSphere>();
	}
	return std::shared_ptr<Scenario>();
}

std::shared_ptr<Scenario> FactorStatic(const std::string& fn)
{
	return std::make_shared<ScenarioStatic>(fn);
}

std::shared_ptr<Scenario> FactorOni(const std::string& fn, unsigned int offset)
{
	return std::make_shared<ScenarioOni>(fn, offset);
}

std::shared_ptr<Scenario> FactorLive()
{
	return std::shared_ptr<Scenario>();
}

