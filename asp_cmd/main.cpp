#include <asp/asp.hpp>
#include <density/PointDensity.hpp>
#include <slimage/image.hpp>
#include <slimage/opencv.hpp>
#include <slimage/io.hpp>
#include <boost/program_options.hpp>
#include <boost/timer/timer.hpp>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <cmath>

float TestDensityFunction(float x, float y)
{
	// constexpr float PI = 3.1415f;
	// float u = 2.0f*x - 1.0f;
	// float v = 2.0f*y - 1.0f;
	// float h = ((1.0f-y)*std::cos(1.5f*PI*u) + y*std::cos(2.5f*PI*u))*std::cos(1.5f*PI*v);
	// return 0.1f + 4.0f*(1.0f+h)*(1.0f+h);
	return 1.0f;
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

Eigen::MatrixXf ComputePointDensity(const asp::Superpixels<Eigen::Vector3f>& sp, const Eigen::MatrixXf& ref)
{
	std::vector<Eigen::Vector2f> seeds(sp.clusters.size());
	std::transform(sp.clusters.begin(), sp.clusters.end(), seeds.begin(),
		[](const asp::Cluster<Eigen::Vector3f>& c) {
			return Eigen::Vector2f{ c.x, c.y };
		});
	return density::PointDensity(seeds, ref);
}

int main(int argc, char** argv)
{
	std::string p_feature_img = "";
	std::string p_density = "";
	std::string p_fn_points = "";
	std::string p_out = "";
	unsigned p_num = 250;
	bool p_verbose = false;
	unsigned p_repetitions = 1;
	bool p_error = false;
	bool p_save_density = true;

	// parameters
	asp::Parameters params = asp::parameters_default();

	namespace po = boost::program_options;
	po::options_description desc;
	desc.add_options()
		("help", "produce help message")
		("features", po::value(&p_feature_img), "feature image (leave empty for uniform image)")
		("density", po::value(&p_density), "density function image (leave empty for test function)")
		("points", po::value(&p_fn_points), "file with points to use as seeds (for DDS)")
		("out", po::value(&p_out), "filename of result file with samples points")
		("num", po::value(&p_num), "number of points to sample")
		("p_num_iterations", po::value(&params.num_iterations), "number of DALIC iterations")
		("p_pds_mode", po::value(&params.pds_mode), "Poisson Disk Sampling method")
		("p_seed_mean_radius_factor", po::value(&params.seed_mean_radius_factor), "size factor for initial cluster mean feature")
		("p_coverage", po::value(&params.coverage), "DALIC cluster search factor")
		("p_weight_compact", po::value(&params.weight_compact), "weight for compactness term")
		("verbose", po::value(&p_verbose), "verbose")
		("repetitions", po::value(&p_repetitions), "repetitions")
		("error", po::value(&p_error), "error")
		("save_density", po::value(&p_save_density), "save_density")
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
		if(p_verbose) std::cout << "Loaded features dim=" << width << "x" << height << "." << std::endl;
	}
	if(!is_density_null) {
		rho = density::LoadDensity(p_density);
		if(p_verbose) std::cout << "Loaded density dim=" << rho.rows() << "x" << rho.cols() << ", sum=" << rho.sum() << "." << std::endl;
	}
	if(is_features_null) {
		width = rho.rows();
		height = rho.cols();
		features.resize(width*height, Eigen::Vector3f{1,1,1});
		if(p_verbose) std::cout << "Created uniform features dim=" << rho.rows() << "x" << rho.cols() << "." << std::endl;
	}
	if(is_density_null) {
		//rho = CreateDensity(p_size, p_size, [](float x, float y) { return TestDensityFunction(x,y); } );
		rho = CreateDensity(width, height, [](float x, float y) { return 1.0f; } );
		if(p_verbose) std::cout << "Created density dim=" << rho.rows() << "x" << rho.cols() << ", sum=" << rho.sum() << "." << std::endl;
	}
	assert(width == rho.rows());
	assert(height == rho.cols());

	std::vector<Eigen::Vector2f> seed_points;
	if(!p_fn_points.empty()) {
		std::ifstream ifs(p_fn_points);
		if(!ifs.is_open()) {
			std::cerr << "Error opening file '" << p_fn_points << "'" << std::endl;
		}
		std::string line;
		while(getline(ifs, line)) {
			std::istringstream ss(line);
			std::vector<float> q;
			while(!ss.eof()) {
				float v;
				ss >> v;
				q.push_back(v);
			}
			seed_points.push_back({q[0], q[1]});
		}
	}

	// scale density
	float scl = static_cast<float>(p_num) / rho.sum();
	rho *= scl;


	// asp
	asp::Superpixels<Eigen::Vector3f> superpixels;
	double total_time = 0.0;
	double total_error = 0.0;
	std::vector<asp::Superpixels<Eigen::Vector3f>> superpixels_v(p_repetitions);
	{
		boost::timer::cpu_timer t;
		for(int k=0; k<p_repetitions; k++) {
			superpixels_v[k] = asp::AdaptiveSuperpixelsRGB(rho, seed_points, features, params);
		}
		auto dt = t.elapsed();
		total_time = static_cast<double>(dt.user + dt.system) / 1000000000.0;
		total_time /= static_cast<double>(p_repetitions);
	}
	if(p_error) {
		for(int k=0; k<p_repetitions; k++) {
			Eigen::MatrixXf approx = ComputePointDensity(superpixels_v[k], rho);
			float error = 0.0;
			for(int i=0; i<approx.size(); i++) {
				error += std::abs(approx.data()[i] - rho.data()[i]);
			}
			error /= rho.sum();
			total_error += error;
		}
		total_error /= static_cast<double>(p_repetitions);
	}
	superpixels = superpixels_v.front();
	if(p_verbose) std::cout << "Generated " << superpixels.clusters.size() << " superpixels." << std::endl;
	std::cout << "T=" << total_time << std::endl;
	std::cout << "E=" << total_error << std::endl;

	// output
	if(!p_out.empty()) {
		// clusters
		{
			std::string fn_clusters = p_out + "_clusters.tsv";
			std::ofstream ofs(fn_clusters);
			for(const auto& c : superpixels.clusters) {
				ofs << c.x << "\t" << c.y << "\t" << c.f[0] << "\t" << c.f[1] << "\t" << c.f[2] << std::endl;
			}
			if(p_verbose) std::cout << "Wrote clusters to file '" << fn_clusters << "'." << std::endl;
		}

		// labels
		if(!superpixels.labels.empty()) {
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
			if(p_verbose) std::cout << "Wrote labels to file '" << fn_labels << "'." << std::endl;
		}

		// commanded density
		if(p_save_density) {
			std::string fn_actual_density = p_out + "_actual_density.tsv";
			density::SaveDensity(fn_actual_density, rho);
			if(p_verbose) std::cout << "Wrote actual density to file '" << fn_actual_density << "'." << std::endl;
		}

		// resulting density
		if(p_save_density) {
			Eigen::MatrixXf approx = ComputePointDensity(superpixels, rho);
			std::string fn_point_density = p_out + "_point_density.tsv";
			density::SaveDensity(fn_point_density, approx);
			if(p_verbose) std::cout << "Wrote point density to file '" << fn_point_density << "'." << std::endl;
		}
	}

	return 1;
}
