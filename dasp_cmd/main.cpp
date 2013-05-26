#include <dasp/Superpixels.hpp>
#include <dasp/Plots.hpp>
#include <dasp/Segmentation.hpp>
#include <Slimage/IO.hpp>
#include <Slimage/Slimage.hpp>
#include <boost/program_options.hpp>
#include <iostream>

int main(int argc, char** argv)
{
	std::string p_img = "";
	std::string p_out = "out";
	bool p_verbose = false;

	dasp::Parameters opt;
	opt.camera = dasp::Camera{320.0f, 240.0f, 540.0f, 0.001f};
	opt.weight_spatial = 1.0f;
	opt.weight_color = 2.0f;
	opt.weight_normal = 3.0f;
	opt.base_radius = 0.018f;
	opt.count = 0;

	namespace po = boost::program_options;
	po::options_description desc;
	desc.add_options()
		("help", "produce help message")
		("img", po::value(&p_img), "path to input image (must have X_color.png and X_depth.pgm)")
		("out", po::value(&p_out), "path to result (will write X.png, X_clusters.tsv and X_labels.tsv)")
		("verbose", po::value(&p_verbose), "verbose")
		("p_radius", po::value(&opt.base_radius), "superpixel radius (meters)")
		("p_count", po::value(&opt.count), "number of superpixels (set to 0 to use radius)")
		("p_num_iterations", po::value(&opt.iterations), "number of DALIC iterations")
		("p_weight_spatial", po::value(&opt.weight_spatial), "weight spatial")
		("p_weight_color", po::value(&opt.weight_color), "weight color")
		("p_weight_normal", po::value(&opt.weight_normal), "weight normal")
	;

	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	po::notify(vm);
	if(vm.count("help")) {
		std::cerr << desc << std::endl;
		return 1;
	}

	const std::string p_img_color = p_img + "_color.png";
	const std::string p_img_depth = p_img + "_depth.pgm";
	if(p_verbose) std::cout << "Reading COLOR: '" << p_img_color << "'" << std::endl;
	if(p_verbose) std::cout << "Reading DEPTH: '" << p_img_depth << "'" << std::endl;
	slimage::Image3ub img_color = slimage::Load3ub(p_img_color);
	slimage::Image1ui16 img_depth =  slimage::Load1ui16(p_img_depth);

	if(p_verbose) std::cout << "n_given=" << opt.count << std::flush;
	dasp::Superpixels superpixels = dasp::ComputeSuperpixels(img_color, img_depth, opt);
	if(p_verbose) std::cout << " R=" << superpixels.opt.base_radius << " n_seeds=" << superpixels.seeds.size() << " n_final=" << superpixels.cluster.size() << std::endl;

	if(!p_out.empty()) {
		const slimage::Image1i& labels = superpixels.ComputeLabels();

		// image
		{
			std::string fn_img = p_out + ".png";
			slimage::Image3ub result = dasp::plots::PlotPoints(superpixels, dasp::plots::Color);
			dasp::plots::PlotEdges(result, labels, slimage::Pixel3ub{{255,0,0}}, 1);
			if(p_verbose) std::cout << "Writing superpixel image to '" << fn_img << "'" << std::endl;
			slimage::Save(result, fn_img);
		}

		// clusters
		{
			std::string fn_clusters = p_out + "_clusters.tsv";
			std::ofstream ofs(fn_clusters);
			for(const auto& cluster : superpixels.cluster) {
				const auto& c = cluster.center;
				ofs
				<< c.px << "\t" << c.py << "\t"
				<< c.color[0] << "\t" << c.color[1] << "\t" << c.color[2] << "\t"
				<< c.position[0] << "\t" << c.position[1] << "\t" << c.position[2] << "\t"
				<< c.normal[0] << "\t" << c.normal[1] << "\t" << c.normal[2] << std::endl;
			}
			if(p_verbose) std::cout << "Wrote clusters to file '" << fn_clusters << "'." << std::endl;
		}

		// labels
		{
			std::string fn_labels = p_out + "_labels.tsv";
			std::ofstream ofs(fn_labels);
			int h = labels.height();
			int w = labels.width();
			for(int y=0; y<h; y++) {
				for(int x=0; x<w; x++) {
					ofs << labels[x + y*w];
					if(x+1 != w) {
						ofs << "\t";
					}
					else {
						ofs << std::endl;
					}
				}
			}
			if(p_verbose) std::cout << "Wrote labels to file '" << fn_labels << "'." << std::endl;
		}
	}

	return 0;
}
