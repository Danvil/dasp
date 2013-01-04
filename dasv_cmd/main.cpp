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

	std::string p_rgbd_mode = "test";
	std::string p_rgbd_arg = "uniform";
	unsigned int p_rgbd_seek = 0;

	namespace po = boost::program_options;
	po::options_description desc;
	desc.add_options()
		("help", "produce help message")
		("rgbd_mode", po::value(&p_rgbd_mode)->default_value(p_rgbd_mode), "rgbd stream mode (test, static, oni, live)")
		("rgbd_arg", po::value(&p_rgbd_arg)->default_value(p_rgbd_arg), "rgbd stream argument")
		("rgbd_seek", po::value(&p_rgbd_seek)->default_value(p_rgbd_seek), "rgbd stream seek (only usable with mode=oni)")
	;

	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	po::notify(vm);
	if(vm.count("help")) {
		std::cerr << desc << std::endl;
		return 1;
	}

	std::cout << "Opening RGBD stream..." << std::endl;
	std::shared_ptr<RgbdStream> rgbd_stream = FactorStream(p_rgbd_mode, p_rgbd_arg);
	auto rgbd_stream_ra = std::dynamic_pointer_cast<RandomAccessRgbdStream>(rgbd_stream);
	if(rgbd_stream_ra) {
		rgbd_stream_ra->seek(p_rgbd_seek);
	}

	std::cout << "Running depth-adaptive supervoxel streaming..." << std::endl;
	dasv::ContinuousSupervoxels sv;
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

	std::cout << "Final supervoxel count = " << sv.numClusters() << std::endl;

	std::cout << "Finished." << std::endl;

	return 1;
}
