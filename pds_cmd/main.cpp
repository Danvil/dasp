#include <pds/PDS.hpp>
#include <density/PointDensity.hpp>
#include <Slimage/Slimage.hpp>
#define SLIMAGE_IO_OPENCV
#include <Slimage/IO.hpp>
#include <boost/program_options.hpp>
#include <boost/timer/timer.hpp>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <cmath>

float TestFunction(float x, float y)
{
	constexpr float PI = 3.1415f;
	float u = 2.0f*x - 1.0f;
	float v = 2.0f*y - 1.0f;
	float h = ((1.0f-y)*std::cos(1.5f*PI*u) + y*std::cos(2.5f*PI*u))*std::cos(1.5f*PI*v);
	return 0.1f + 4.0f*(1.0f+h)*(1.0f+h);
}

Eigen::MatrixXf TestDensity(unsigned size, unsigned num)
{
	Eigen::MatrixXf rho(size,size);
	for(unsigned i=0; i<size; i++) {
		float y = static_cast<float>(i) / static_cast<float>(size-1);
		for(unsigned j=0; j<size; j++) {
			float x = static_cast<float>(j) / static_cast<float>(size-1);
			rho(j,i) = TestFunction(x,y);
		}
	}
	return rho;
}

int main(int argc, char** argv)
{
	std::string p_density = "";
	std::string p_mode = "spds";
	std::string p_out = "pnts.tsv";
	unsigned p_size = 128;
	unsigned p_num = 250;

	namespace po = boost::program_options;
	po::options_description desc;
	desc.add_options()
		("help", "produce help message")
		("density", po::value(&p_density), "density function image (leave empty for test function)")
		("mode", po::value(&p_mode), "sampling method")
		("out", po::value(&p_out), "filename of result file with samples points")
		("size", po::value(&p_size), "size of image in pixel")
		("num", po::value(&p_num), "number of points to sample")
	;

	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	po::notify(vm);
	if(vm.count("help")) {
		std::cerr << desc << std::endl;
		return 1;
	}

	Eigen::MatrixXf rho;
	if(p_density.empty()) {
		rho = TestDensity(p_size, p_num);
		std::cout << "Created density dim=" << rho.rows() << "x" << rho.cols() << ", sum=" << rho.sum() << "." << std::endl;
	}
	else {
		rho = density::LoadDensity(p_density);
		std::cout << "Loaded density dim=" << rho.rows() << "x" << rho.cols() << ", sum=" << rho.sum() << "." << std::endl;
	}

	// scale density
	float scl = static_cast<float>(p_num) / rho.sum();
	rho *= scl;

	std::vector<Eigen::Vector2f> pnts;
	{
		boost::timer::auto_cpu_timer t;
		pnts = pds::PoissonDiscSampling(p_mode, rho);
	}
	std::cout << "Generated " << pnts.size() << " points." << std::endl;

	std::ofstream ofs(p_out);
	for(const Eigen::Vector2f& x : pnts) {
		ofs << x[0] << "\t" << x[1] << std::endl;
	}
	std::cout << "Wrote points to file '" << p_out << "'." << std::endl;

	return 1;
}
