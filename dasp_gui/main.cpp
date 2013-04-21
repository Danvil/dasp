#include "WdgtKinectSuperPoints.h"
#include <QtGui>
#include <QApplication>
#include <boost/program_options.hpp>
#include <iostream>
#include <string>

int main(int argc, char *argv[])
{
	namespace po = boost::program_options;
	bool p_no3d = false;

	// parse command line options
	po::options_description desc("Allowed options");
	desc.add_options()
		("help", "produce help message")
		("no3d", po::value(&p_no3d), "Disable 3D rendering")
		("live", "Kinect live mode")
		("oni", po::value<std::string>(), "view a pre-recorded ONI file")
		("rgbd", po::value<std::string>(), "open an RGB-D image")
	;

	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	po::notify(vm);    

	if (vm.count("help")) {
		std::cout << desc << std::endl;
		return 1;
	}

    QApplication a(argc, argv);
    WdgtKinectSuperPoints w(p_no3d);
    w.show();

	if(vm.count("live")) {
		w.ShowLive();
	}
	else if(vm.count("oni")) {
		w.LoadOni(vm["oni"].as<std::string>());
	}
	else if(vm.count("rgbd")) {
		w.LoadRgbd(vm["rgbd"].as<std::string>());
	}

    return a.exec();
}
