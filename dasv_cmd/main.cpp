#include <dasv.hpp>
#include <Slimage/IO.hpp>
#include <Slimage/Gui.hpp>
#include <kinect/KinectGrabber.h>
#include <boost/program_options.hpp>
#include <string>
#include <iostream>
#include <memory>

constexpr int WIDTH = 640;
constexpr int HEIGHT = 480;

struct Data
{
	slimage::Image3ub color;
	slimage::Image1ui16 depth;
};

class Scenario
{
public:
	virtual ~Scenario() {}
	virtual Data pop() = 0;
};

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
	Data pop() { return data_; }
private:
	Data data_;
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
	Data pop() { return data_; }
private:
	Data data_;
};

class ScenarioTestSphere : public Scenario
{
public:
	ScenarioTestSphere() : time_(0) {}
	Data pop() {
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
	Data pop() { return data_; }
private:
	Data data_;
};

class ScenarioOni : public Scenario
{
public:
	ScenarioOni(const std::string& fn, int offset=0) {
		grabber_ = std::make_shared<dasp::KinectGrabber>();
		grabber_->OpenFile(fn);
		grabber_->SeekToFrame(offset);
	}
	Data pop() {
		Data data;
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

int main(int argc, char** argv)
{
	const std::string ds_path = "/home/david/Documents/DataSets";

	std::string p_mode;

	namespace po = boost::program_options;
	po::options_description desc;
	desc.add_options()
		("help", "produce help message")
		("mode", po::value<std::string>(&p_mode)->required(), "mode string")
	;

	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	po::notify(vm);
	if(vm.count("help")) {
		std::cerr << desc << std::endl;
		return 1;
	}

	using namespace dasv;

	std::shared_ptr<Scenario> scn;

	if(p_mode == "test_uniform") {
		scn.reset(new ScenarioTestUniform());
	}
	else if(p_mode == "test_paraboloid") {
		scn.reset(new ScenarioTestParaboloid());
	}
	else if(p_mode == "test_sphere") {
		scn.reset(new ScenarioTestSphere());
	}
	else if(p_mode == "static") {
		std::string fn = ds_path + "/dasp_rgbd_dataset/images/001";
		scn.reset(new ScenarioStatic(fn));
	}
	else if(p_mode == "oni") {
		std::string fn = ds_path + "/2012-10-12 cogwatch dasp/SlowVelocity/C15_c01_slow.oni";
		scn.reset(new ScenarioOni(fn, 100));
	}

	ContinuousSupervoxels sv;
	sv.start(WIDTH,HEIGHT);
	while(true) {
		Data data = scn->pop();
		if(!data.color || !data.depth) {
			break;
		}
		slimage::gui::Show("color", data.color, 0);
		slimage::gui::Show("depth", data.depth, 500, 3000, 0);
		sv.step(data.color, data.depth);
	}

	std::cout << "Supervoxel count = " << sv.numClusters() << std::endl;

	std::cout << "Finished." << std::endl;

	return 1;
}
