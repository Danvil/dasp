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

bool p_verbose = false;

namespace impl
{
	void write_result(std::ostream& ofs, const std::vector<float>& v)
	{
		for(unsigned int i=0; i<v.size(); i++) {
			if(i > 0)
				ofs << "," << std::endl;
			ofs << v[i] << std::endl;
		}
	}

	template<typename F>
	std::vector<std::vector<float>> evaluate(
		const slimage::Image3ub& img_color, const slimage::Image1ui16& img_depth,
		const dasp::Parameters& opt,
		unsigned int num,
		F f
	)
	{
		std::vector<std::vector<float>> result;
		float total = 0.0f;
		for(unsigned int k=0; k<num; k++) {
			dasp::Superpixels dasp = dasp::ComputeSuperpixels(img_color, img_depth, opt);
			auto v = f(dasp);
			result.push_back(v);
		}
		return result;
	}

	std::vector<float> mean(const std::vector<std::vector<float>>& v) {
		std::vector<float> u(v[0].size(), 0.0f);
		for(unsigned int i=0; i<v.size(); i++) {
			assert(v[i].size() == u.size());
			for(unsigned int j=0; j<u.size(); j++) {
				u[j] += v[i][j];
			}
		}
		for(unsigned int j=0; j<u.size(); j++) {
			u[j] /= static_cast<float>(v.size());
		}
		return u;
	}
}

template<typename F>
void process(
	const std::string& mode, const std::string& name, const std::string& filename,
	const slimage::Image3ub& img_color, const slimage::Image1ui16& img_depth,
	const dasp::Parameters& opt,
	unsigned int num,
	F f
)
{
	auto q = impl::mean(impl::evaluate(img_color, img_depth, opt, num, f));
	if(p_verbose || filename.empty()) {
		std::cout << name << ": ";
		impl::write_result(std::cout, q);
	}
	if(!filename.empty()) {
		std::ofstream ofs(filename);
		ofs << mode << ",";
		impl::write_result(ofs, q);
	}
}

dasp::DensityMode StringToDensityMode(const std::string& dm)
{
	if(dm == "ASP_const") return dasp::DensityMode::ASP_const;
	if(dm == "ASP_depth") return dasp::DensityMode::ASP_depth;
	if(dm == "DASP") return dasp::DensityMode::DASP;
	exit(0);
}

int main(int argc, char** argv)
{
	std::string p_mode;
	std::string p_density = "DASP";
	std::string p_img_path;
	std::string p_truth_path;
	std::string p_result_path;
	unsigned int p_br_d = 2;
	unsigned int p_num = 5;

	dasp::Parameters opt;
	opt.camera = dasp::Camera{320.0f, 240.0f, 540.0f, 0.001f};
	opt.weight_spatial = 1.0f;
	opt.weight_color = 2.0f;
	opt.weight_normal = 3.0f;

	po::options_description desc;
	desc.add_options()
		("help", "produce help message")
		("mode", po::value<std::string>(&p_mode), "operating mode: use, br, ...")
		("density_mode", po::value<std::string>(&p_density), "density mode: ASP_const, ASP_depth, DASP")
		("image", po::value<std::string>(&p_img_path), "path to RGBD image")
		("truth", po::value<std::string>(&p_truth_path), "path to ground truth")
		("result", po::value<std::string>(&p_result_path), "path to result")
		("radius", po::value(&opt.base_radius)->default_value(opt.base_radius), "superpixel radius (meters)")
		("count", po::value(&opt.count)->default_value(opt.count), "number of superpixels (set to 0 to use radius)")
		("iterations", po::value(&opt.iterations)->default_value(opt.iterations), "number of iterations for local nearest neighbour clustering")
		("repetitions", po::value(&p_num)->default_value(p_num), "number of repetitions")
		("br_d", po::value(&p_br_d)->default_value(p_br_d), "border distance tolerance in pixel")
	;

	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	po::notify(vm);
	if(vm.count("help")) {
		std::cerr << desc << std::endl;
		return 1;
	}

	opt.density_mode = StringToDensityMode(p_density);

	const std::string p_img_path_color = p_img_path + "_color.png";
	const std::string p_img_path_depth = p_img_path + "_depth.pgm";
	if(p_verbose) std::cout << "Reading COLOR: '" << p_img_path_color << "'" << std::endl;
	if(p_verbose) std::cout << "Reading DEPTH: '" << p_img_path_depth << "'" << std::endl;
	slimage::Image3ub img_color = slimage::Load3ub(p_img_path_color);
	slimage::Image1ui16 img_depth =  slimage::Load1ui16(p_img_path_depth);

	if(p_mode == "use") {
		// ground truth labels
		if(p_verbose) std::cout << "Reading TRUTH: '" << p_truth_path << "'" << std::endl;
		slimage::Image1i img_truth;
		slimage::conversion::Convert(slimage::Load1ui16(p_truth_path), img_truth);
		// eval
		process("use", "Undersegmentation error (USE)", p_result_path,
			img_color, img_depth, opt, p_num,
			[=](const dasp::Superpixels& superpixels) -> std::vector<float> {
				// actual labels
				slimage::Image1i labels = superpixels.ComputeLabels();
				// undersegmentation error
				float use = dasp::eval::UndersegmentationError(img_truth, labels);
				return { use };
			});
	}

	if(p_mode == "br") {
		// ground truth boundaries
		if(p_verbose) std::cout << "Reading TRUTH: '" << p_truth_path << "'" << std::endl;
		slimage::Image1ub truth_boundary = slimage::Pick<unsigned char>(slimage::Load(p_truth_path), 0);
		// eval
		process("br", "Boundary recall (BR)", p_result_path,
			img_color, img_depth, opt, p_num,
			[=](const dasp::Superpixels& superpixels) -> std::vector<float> {
				// actual boundaries
				slimage::Image1i labels = superpixels.ComputeLabels();
				slimage::Image1ub img_boundaries(superpixels.width(), superpixels.height(), slimage::Pixel1ub{0});
				dasp::plots::PlotEdges(img_boundaries, labels, slimage::Pixel1ub{255}, 2);
				// boundary recall
				return {dasp::ComputeRecallBox(truth_boundary, img_boundaries, p_br_d)};
			});
	}

	if(p_mode == "area") {
		process("area", "Superpixel Area (AREA)", p_result_path,
			img_color, img_depth, opt, p_num,
			[](const dasp::Superpixels& superpixel) -> std::vector<float> {
				return { dasp::eval::Area(superpixel) };
			});
	}

	if(p_mode == "area3") {
		process("area3", "Superpixel Area 3D (AREA3)", p_result_path,
			img_color, img_depth, opt, p_num,
			[](const dasp::Superpixels& superpixel) -> std::vector<float> {
				return { dasp::eval::Area3D(superpixel) };
			});
	}

	if(p_mode == "ipq") {
		process("ipq", "Isoperimetric Quotient (IPQ)", p_result_path,
			img_color, img_depth, opt, p_num,
			[](const dasp::Superpixels& superpixel) -> std::vector<float> {
				std::pair<float,std::vector<float>> ipq = dasp::eval::IsoperimetricQuotient(superpixel);
				return { ipq.first };
			});
	}

	if(p_mode == "ipq3") {
		process("ipq3", "Isoperimetric Quotient 3D (IPQ3)", p_result_path,
			img_color, img_depth, opt, p_num,
			[](const dasp::Superpixels& superpixel) -> std::vector<float> {
				std::pair<float,std::vector<float>> ipq = dasp::eval::IsoperimetricQuotient3D(superpixel);
				return { ipq.first };
			});
	}

	if(p_mode == "ev_c") {
		process("ev_c", "Explained Variation (EV) - Color", p_result_path,
			img_color, img_depth, opt, p_num,
			[](const dasp::Superpixels& superpixel) -> std::vector<float> {
				return { dasp::eval::ExplainedVariationColor(superpixel) };
			});
	}

	if(p_mode == "ev_d") {
		process("ev_d", "Explained Variation (EV) - Depth", p_result_path,
			img_color, img_depth, opt, p_num,
			[](const dasp::Superpixels& superpixel) -> std::vector<float> {
				return { dasp::eval::ExplainedVariationDepth(superpixel) };
			});
	}

	if(p_mode == "ev_v") {
		process("ev_v", "Explained Variation (EV) - Position", p_result_path,
			img_color, img_depth, opt, p_num,
			[](const dasp::Superpixels& superpixel) -> std::vector<float> {
				return { dasp::eval::ExplainedVariationPosition(superpixel) };
			});
	}

	if(p_mode == "ev_n") {
		process("ev_n", "Explained Variation (EV) - Normal", p_result_path,
			img_color, img_depth, opt, p_num,
			[](const dasp::Superpixels& superpixel) -> std::vector<float> {
				return { dasp::eval::ExplainedVariationNormal(superpixel) };
			});
	}

	if(p_mode == "ce_c") {
		process("ce_c", "Compression Error (CE) - Color", p_result_path,
			img_color, img_depth, opt, p_num,
			[](const dasp::Superpixels& superpixel) -> std::vector<float> {
				return { dasp::eval::CompressionErrorColor(superpixel) };
			});
	}

	if(p_mode == "ce_d") {
		process("ce_d", "Compression Error (CE) - Depth", p_result_path,
			img_color, img_depth, opt, p_num,
			[](const dasp::Superpixels& superpixel) -> std::vector<float> {
				return { dasp::eval::CompressionErrorDepth(superpixel) };
			});
	}

	if(p_mode == "ce_v") {
		process("ce_v", "Compression Error (CE) - Position", p_result_path,
			img_color, img_depth, opt, p_num,
			[](const dasp::Superpixels& superpixel) -> std::vector<float> {
				return { dasp::eval::CompressionErrorPosition(superpixel) };
			});
	}

	if(p_mode == "ce_n") {
		process("ce_n", "Compression Error (CE) - Normal", p_result_path,
			img_color, img_depth, opt, p_num,
			[](const dasp::Superpixels& superpixel) -> std::vector<float> {
				return { dasp::eval::CompressionErrorNormal(superpixel) };
			});
	}

	if(p_mode == "nd") {
		process("nd", "Mean Neighbour Distance (MND)", p_result_path,
			img_color, img_depth, opt, p_num,
			[](const dasp::Superpixels& superpixel) -> std::vector<float> {
				return { dasp::eval::MeanNeighbourDistance(superpixel) };
			});
	}

	if(p_mode == "ew_thick") {
		process("ew_thick", "Eigenvalues Thickness", p_result_path,
			img_color, img_depth, opt, p_num,
			[](const dasp::Superpixels& superpixel) -> std::vector<float> {
				float q = std::accumulate(superpixel.cluster.begin(), superpixel.cluster.end(), 0.0f,
					[](float a, const dasp::Cluster& c) {
						return a + c.thickness;
					}) / static_cast<float>(superpixel.cluster.size());
				return { q };
			});
	}

	if(p_mode == "ew_ecc") {
		process("ew_ecc", "Eigenvalues Eccentricity", p_result_path,
			img_color, img_depth, opt, p_num,
			[](const dasp::Superpixels& superpixel) -> std::vector<float> {
				// for(unsigned int i=0; i<superpixel.cluster.size(); i++) {
				// 	std::cout << superpixel.cluster[i].eccentricity << std::endl;
				// }
				float q = std::accumulate(superpixel.cluster.begin(), superpixel.cluster.end(), 0.0f,
					[](float a, const dasp::Cluster& c) {
						return a + c.eccentricity;
					}) / static_cast<float>(superpixel.cluster.size());
				return { q };
			});
	}

	if(p_mode == "ew_aq") {
		process("ew_aq", "Eigenvalues Area Quotient", p_result_path,
			img_color, img_depth, opt, p_num,
			[](const dasp::Superpixels& superpixel) -> std::vector<float> {
				float q = std::accumulate(superpixel.cluster.begin(), superpixel.cluster.end(), 0.0f,
					[](float a, const dasp::Cluster& c) {
						return a + c.area_quotient;
					}) / static_cast<float>(superpixel.cluster.size());
				return { q };
			});
	}

	if(p_mode == "ew_area") {
		process("ew_area", "Eigenvalues Area", p_result_path,
			img_color, img_depth, opt, p_num,
			[](const dasp::Superpixels& superpixel) -> std::vector<float> {
				float q = std::accumulate(superpixel.cluster.begin(), superpixel.cluster.end(), 0.0f,
					[](float a, const dasp::Cluster& c) {
						return a + c.area;
					}) / static_cast<float>(superpixel.cluster.size());
				return { q };
			});
	}

	return 0;
}
