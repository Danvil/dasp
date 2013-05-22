#include <asp/asp.hpp>
#include <Slimage/Slimage.hpp>
#define SLIMAGE_IO_OPENCV
#include <Slimage/IO.hpp>
#include <boost/program_options.hpp>
#include <boost/timer/timer.hpp>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <cmath>

float TestDensityFunction(float x, float y)
{
	constexpr float PI = 3.1415f;
	float u = 2.0f*x - 1.0f;
	float v = 2.0f*y - 1.0f;
	float h = ((1.0f-y)*std::cos(1.5f*PI*u) + y*std::cos(2.5f*PI*u))*std::cos(1.5f*PI*v);
	return 0.1f + 4.0f*(1.0f+h)*(1.0f+h);
}

template<typename F>
Eigen::MatrixXf CreateDensity(unsigned width, unsigned height, F f)
{
	Eigen::MatrixXf rho(width,height);
	for(unsigned i=0; i<height; i++) {
		float y = static_cast<float>(i) / static_cast<float>(height-1);
		for(unsigned j=0; j<width; j++) {
			float x = static_cast<float>(j) / static_cast<float>(width-1);
			rho(j,i) = f(x,y);
		}
	}
	return rho;
}

Eigen::MatrixXf LoadDensity(const std::string& filename)
{
	slimage::ImagePtr ptr = slimage::Load(filename);
	if(!ptr) {
		std::cerr << "Could not load image!" << std::endl;
		throw 0;
	}
	if(slimage::HasType<unsigned char, 3>(ptr)) {
		slimage::Image3ub img = slimage::Ref<unsigned char, 3>(ptr);
		Eigen::MatrixXf mat(img.width(), img.height());
		for(int y=0; y<mat.cols(); y++) {
			for(int x=0; x<mat.rows(); x++) {
				slimage::Pixel3ub p = img(x,y);
				mat(x,y) = static_cast<float>((int)p[0] + (int)p[1] + (int)p[2])/3.0f/255.0f;
			}
		}
		return mat;
	}
	if(slimage::HasType<unsigned char, 1>(ptr)) {
		slimage::Image1ub img = slimage::Ref<unsigned char, 1>(ptr);
		Eigen::MatrixXf mat(img.width(), img.height());
		for(int y=0; y<mat.cols(); y++) {
			for(int x=0; x<mat.rows(); x++) {
				mat(x,y) = static_cast<float>(img(x,y))/255.0f;
			}
		}
		return mat;
	}
}

std::vector<Eigen::Vector3f> LoadFeatures(const std::string& fn, int& width, int& height)
{
	slimage::Image3ub img = slimage::Load3ub(fn);
	std::vector<Eigen::Vector3f> features;
	features.resize(img.width()*img.height());
	for(int y=0; y<img.height(); y++) {
		for(int x=0; x<img.width(); x++) {
			slimage::Pixel3ub p = img(x,y);
			Eigen::Vector3f col(
				static_cast<float>(p[0])/255.0f,
				static_cast<float>(p[1])/255.0f,
				static_cast<float>(p[2])/255.0f
				);
			features[x + y*img.width()] = col;
		}
	}
	width = img.width();
	height = img.height();
	return features;
}

int main(int argc, char** argv)
{
	std::string p_feature_img = "";
	std::string p_density = "";
	std::string p_out = "sp";
	unsigned p_num = 250;

	namespace po = boost::program_options;
	po::options_description desc;
	desc.add_options()
		("help", "produce help message")
		("features", po::value(&p_feature_img), "feature image (leave empty for uniform image)")
		("density", po::value(&p_density), "density function image (leave empty for test function)")
		("out", po::value(&p_out), "filename of result file with samples points")
		("num", po::value(&p_num), "number of points to sample")
	;

	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	po::notify(vm);
	if(vm.count("help")) {
		std::cerr << desc << std::endl;
		return 1;
	}


	// load

	bool is_features_null = p_feature_img.empty();
	bool is_density_null = p_density.empty();
	
	Eigen::MatrixXf rho;
	std::vector<Eigen::Vector3f> features;
	int width = -1, height = -1;

	if(is_density_null && is_features_null) {
		std::cerr << "ERROR feature or density must not be both empty!" << std::endl;
		return 1;
	}
	if(!is_features_null) {
		features = LoadFeatures(p_feature_img, width, height);
		std::cout << "Loaded features dim=" << width << "x" << height << "." << std::endl;
	}
	if(!is_density_null) {
		rho = LoadDensity(p_density);
		std::cout << "Loaded density dim=" << rho.rows() << "x" << rho.cols() << ", sum=" << rho.sum() << "." << std::endl;
	}
	if(is_features_null) {
		features.resize(rho.rows()*rho.cols(), Eigen::Vector3f::Zero());
		std::cout << "Created uniform features dim=" << rho.rows() << "x" << rho.cols() << "." << std::endl;
	}
	if(is_density_null) {
		//rho = CreateDensity(p_size, p_size, [](float x, float y) { return TestDensityFunction(x,y); } );
		rho = CreateDensity(width, height, [](float x, float y) { return 1.0f; } );
		std::cout << "Created density dim=" << rho.rows() << "x" << rho.cols() << ", sum=" << rho.sum() << "." << std::endl;
	}
	assert(width == rho.rows());
	assert(height == rho.cols());

	// scale density
	float scl = static_cast<float>(p_num) / rho.sum();
	rho *= scl;


	// parameters
	asp::Parameters params = asp::parameters_default();

	// asp
	asp::Superpixels<Eigen::Vector3f> superpixels = asp::AdaptiveSuperpixelsRGB(
		rho, features, params);
	std::cout << "Generated " << superpixels.clusters.size() << " superpixels." << std::endl;


	// output

	// clusters
	{
		std::string fn_clusters = p_out + "_clusters.tsv";
		std::ofstream ofs(fn_clusters);
		for(const auto& c : superpixels.clusters) {
			ofs << c.x << "\t" << c.y << "\t" << c.f[0] << "\t" << c.f[1] << "\t" << c.f[2] << std::endl;
		}
		std::cout << "Wrote clusters to file '" << fn_clusters << "'." << std::endl;
	}

	// clusters
	{
		std::string fn_labels = p_out + "_labels.tsv";
		std::ofstream ofs(fn_labels);
		for(int y=0; y<height; y++) {
			for(int x=0; x<width; x++) {
				ofs << superpixels.labels[x + y*width];
				if(x+1 != width) {
					ofs << "\t";
				}
				else {
					ofs << std::endl;
				}
			}
		}
		std::cout << "Wrote labels to file '" << fn_labels << "'." << std::endl;
	}

	return 1;
}
