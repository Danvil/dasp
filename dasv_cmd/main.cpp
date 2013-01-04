#include <dasv.hpp>
#include <rgbd/rgbd.hpp>
#include <Slimage/Gui.hpp>
#include <boost/program_options.hpp>
#include <string>
#include <iostream>
#include <memory>

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

	std::shared_ptr<RgbdStream> rgbd_stream;

	if(p_mode == "test_uniform") {
		rgbd_stream = FactorTest("uniform");
	}
	else if(p_mode == "test_paraboloid") {
		rgbd_stream = FactorTest("paraboloid");
	}
	else if(p_mode == "test_sphere") {
		rgbd_stream = FactorTest("sphere");
	}
	else if(p_mode == "static") {
		std::string fn = ds_path + "/dasp_rgbd_dataset/images/001";
		rgbd_stream = FactorStatic(fn);
	}
	else if(p_mode == "oni") {
		std::string fn = ds_path + "/2012-10-12 cogwatch dasp/SlowVelocity/C15_c01_slow.oni";
		auto tmp = FactorOni(fn);
		tmp->seek(100);
		rgbd_stream = tmp;
	}

	ContinuousSupervoxels sv;
	sv.start();
	while(rgbd_stream->grab()) {
		Rgbd data = rgbd_stream->get();
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
