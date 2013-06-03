#include <density/PointDensity.hpp>
#include <density/Smooth.hpp>
#include <boost/program_options.hpp>
#include <boost/timer/timer.hpp>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <cmath>

int main(int argc, char** argv)
{
	std::string p_in = "";
	std::string p_out = "out.tsv";

	namespace po = boost::program_options;
	po::options_description desc;
	desc.add_options()
		("help", "produce help message")
		("density", po::value(&p_in), "filename for input density")
		("out", po::value(&p_out), "filename for smoothed output density")
	;

	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	po::notify(vm);
	if(vm.count("help")) {
		std::cerr << desc << std::endl;
		return 1;
	}

	Eigen::MatrixXf rho = density::LoadDensity(p_in);
	std::cout << "Loaded density dim=" << rho.rows() << "x" << rho.cols() << ", sum=" << rho.sum() << "." << std::endl;

	Eigen::MatrixXf result;
	{
		boost::timer::auto_cpu_timer t;
		result = density::DensityAdaptiveSmooth(rho);
	}

	density::SaveDensity(p_out, result);
	std::cout << "Wrote result to file '" << p_out << "'." << std::endl;

	return 1;
}
