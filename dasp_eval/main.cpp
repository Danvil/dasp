#include <dasp/eval/Recall.hpp>
#include <dasp/Superpixels.hpp>
#include <dasp/Plots.hpp>
#include <dasp/Segmentation.hpp>
#include <dasp/eval/eval.hpp>
#include <Slimage/IO.hpp>
#include <Slimage/Slimage.hpp>
#include <Slimage/Convert.hpp>
#include <boost/program_options.hpp>
namespace po = boost::program_options;
#include <iostream>
#include <fstream>

void write_result(const std::string& filename, const dasp::Superpixels& dasp, const std::string& mode, float v)
{
	std::ofstream ofs(filename);
	ofs << "radius," << dasp.opt.base_radius << std::endl;
	ofs << "num," << dasp.cluster.size() << std::endl;
	ofs << mode << "," << v << std::endl;
}

int main(int argc, char** argv)
{
	std::string p_mode;
	std::string p_img_path;
	std::string p_truth_path;
	std::string p_result_path;
	unsigned int p_br_d = 2;
	bool p_verbose = false;

	dasp::Parameters opt;
	opt.camera = dasp::Camera{320.0f, 240.0f, 540.0f, 0.001f};
	opt.weight_spatial = 1.0f;
	opt.weight_color = 2.0f;
	opt.weight_normal = 3.0f;

	po::options_description desc;
	desc.add_options()
		("help", "produce help message")
		("mode", po::value<std::string>(&p_mode), "operating mode: use, br, ...")
		("image", po::value<std::string>(&p_img_path), "path to RGBD image")
		("truth", po::value<std::string>(&p_truth_path), "path to ground truth")
		("result", po::value<std::string>(&p_result_path)->default_value("/tmp/result.txt"), "path to result")
		("radius", po::value(&opt.base_radius)->default_value(opt.base_radius), "superpixel radius (meters)")
		("count", po::value(&opt.count)->default_value(opt.count), "number of superpixels (set to 0 to use radius)")
		("iterations", po::value(&opt.iterations)->default_value(opt.iterations), "number of iterations for local nearest neighbour clustering")
		("br_d", po::value(&p_br_d)->default_value(p_br_d), "border distance tolerance in pixel")
	;

	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	po::notify(vm);
	if(vm.count("help")) {
		std::cerr << desc << std::endl;
		return 1;
	}

	const std::string p_img_path_color = p_img_path + "_color.png";
	const std::string p_img_path_depth = p_img_path + "_depth.pgm";
	if(p_verbose) std::cout << "Reading COLOR: '" << p_img_path_color << "'" << std::endl;
	if(p_verbose) std::cout << "Reading DEPTH: '" << p_img_path_depth << "'" << std::endl;
	slimage::Image3ub img_color = slimage::Load3ub(p_img_path_color);
	slimage::Image1ui16 img_depth =  slimage::Load1ui16(p_img_path_depth);

	dasp::Superpixels superpixels = dasp::ComputeSuperpixels(img_color, img_depth, opt);
	if(p_verbose) std::cout << " R=" << superpixels.opt.base_radius << " n_seeds=" << superpixels.seeds.size() << " n_final=" << superpixels.cluster.size() << std::endl;
	
	if(p_mode == "use") {
		// ground truth labels
		if(p_verbose) std::cout << "Reading TRUTH: '" << p_truth_path << "'" << std::endl;
		slimage::Image1i img_truth;
		slimage::conversion::Convert(slimage::Load1ui16(p_truth_path), img_truth);
		// actual labels
		slimage::Image1i labels = superpixels.ComputeLabels();
		// undersegmentation error
		std::pair<float,unsigned int> use = dasp::UndersegmentationErrorTotal(img_truth, labels);
		if(p_verbose) std::cout << "Undersegmentation error (USE): " << use.first << " - " << use.second << std::endl;
		write_result(p_result_path, superpixels, p_mode, use.first);
	}

	if(p_mode == "br") {
		// ground truth boundaries
		if(p_verbose) std::cout << "Reading TRUTH: '" << p_truth_path << "'" << std::endl;
		slimage::Image1ub truth_boundary = slimage::Pick<unsigned char>(slimage::Load(p_truth_path), 0);
		// actual boundaries
		slimage::Image1i labels = superpixels.ComputeLabels();
		slimage::Image1ub img_boundaries(superpixels.width(), superpixels.height(), slimage::Pixel1ub{0});
		dasp::plots::PlotEdges(img_boundaries, labels, slimage::Pixel1ub{255}, 2);
		// boundary recall
		float br = dasp::ComputeRecallBox(truth_boundary, img_boundaries, p_br_d);
		if(p_verbose) std::cout << "Boundary recall (BR): " << br << std::endl;
		write_result(p_result_path, superpixels, p_mode, br);
	}

	if(p_mode == "ipq") {
		std::pair<float,std::vector<float>> ipq = dasp::eval::IsoperimetricQuotient(superpixels);
		if(p_verbose) std::cout << "Isoperimetric Quotient (IPQ): " << ipq.first << std::endl;
		write_result(p_result_path, superpixels, p_mode, ipq.first);
	}

	if(p_mode == "ev_c") {
		float ev = dasp::eval::ExplainedVariationColor(superpixels);
		if(p_verbose) std::cout << "Explained Variation (EV) - Color: " << ev << std::endl;
		write_result(p_result_path, superpixels, p_mode, ev);
	}

	if(p_mode == "ev_d") {
		float ev = dasp::eval::ExplainedVariationDepth(superpixels);
		if(p_verbose) std::cout << "Explained Variation (EV) - Depth: " << ev << std::endl;
		write_result(p_result_path, superpixels, p_mode, ev);
	}

	if(p_mode == "ev_v") {
		float ev = dasp::eval::ExplainedVariationPosition(superpixels);
		if(p_verbose) std::cout << "Explained Variation (EV) - Position: " << ev << std::endl;
		write_result(p_result_path, superpixels, p_mode, ev);
	}

	if(p_mode == "ev_n") {
		float ev = dasp::eval::ExplainedVariationNormal(superpixels);
		if(p_verbose) std::cout << "Explained Variation (EV) - Normal: " << ev << std::endl;
		write_result(p_result_path, superpixels, p_mode, ev);
	}

//	const slimage::Image1i& labels = superpixels.ComputeLabels();
//	slimage::Image3ub result = dasp::plots::PlotPoints(superpixels, dasp::plots::Color);
//	dasp::plots::PlotEdges(result, labels, slimage::Pixel3ub{{255,0,0}}, 1);
//	std::cout << "Writing result to '" << p_result_path << "'" << std::endl;
//	slimage::Save(result, p_result_path);

	return 0;
}
