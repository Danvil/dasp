#include "WdgtDasvVis.h"
#include <QtGui>
#include <QApplication>
#include <boost/program_options.hpp>
#include <iostream>
#include <string>

int main(int argc, char *argv[])
{
	namespace po = boost::program_options;

	std::string p_rgbd_mode = "test";
	std::string p_rgbd_arg = "uniform";
	unsigned int p_rgbd_seek = 0;

	// parse command line options
	po::options_description desc("Allowed options");
	desc.add_options()
		("help", "produce help message")
		("rgbd_mode", po::value(&p_rgbd_mode)->default_value(p_rgbd_mode), "rgbd stream mode (test, static, oni, live)")
		("rgbd_arg", po::value(&p_rgbd_arg)->default_value(p_rgbd_arg), "rgbd stream argument")
		("rgbd_seek", po::value(&p_rgbd_seek)->default_value(p_rgbd_seek), "rgbd stream seek (only usable with mode=oni)")
	;

	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	po::notify(vm);    

	if (vm.count("help")) {
		std::cout << desc << std::endl;
		return 1;
	}

	std::cout << "Creating widget..." << std::endl;
	QApplication a(argc, argv);
	WdgtDasvVis w;
	w.show();

	std::cout << "Opening RGBD stream..." << std::endl;
	std::shared_ptr<RgbdStream> rgbd_stream = FactorStream(p_rgbd_mode, p_rgbd_arg);
	auto rgbd_stream_ra = std::dynamic_pointer_cast<RandomAccessRgbdStream>(rgbd_stream);
	if(rgbd_stream_ra) {
		rgbd_stream_ra->seek(p_rgbd_seek);
	}
    w.setRgbdStream(rgbd_stream);

    return a.exec();
}
