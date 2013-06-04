#include <dasp/Superpixels.hpp>
#include <dasp/Plots.hpp>
#include <dasp/Segmentation.hpp>
#include <rgbd/rgbd.hpp>
#include <density/PointDensity.hpp>
#include <Slimage/IO.hpp>
#include <Slimage/Slimage.hpp>
#include <boost/program_options.hpp>
#include <boost/format.hpp>
#include <boost/progress.hpp>
#include <iostream>

int main(int argc, char** argv)
{
	std::string p_img = "";
	std::string p_rgbd_mode = "";
	std::string p_rgbd_arg = "";
	std::string p_out = "out";
	bool p_verbose = false;
	int p_stream_max = 100;
	std::string p_pds_mode = "spds";
	bool p_save_density = false;

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
		("rgbd_mode", po::value(&p_rgbd_mode), "rgbd image stream mode (static, oni, live, ...)")
		("rgbd_arg", po::value(&p_rgbd_arg), "rgbd image stream argument (depends on mode)")
		("out", po::value(&p_out), "path to result (will write X.png, X_clusters.tsv and X_labels.tsv)")
		("verbose", po::value(&p_verbose), "verbose")
		("stream_max", po::value(&p_stream_max), "maximum number of frames to process")
		("p_radius", po::value(&opt.base_radius), "superpixel radius (meters)")
		("p_count", po::value(&opt.count), "number of superpixels (set to 0 to use radius)")
		("p_pds_mode", po::value(&p_pds_mode), "Poisson Disk sampling method (spds, delta)")
		("p_num_iterations", po::value(&opt.iterations), "number of DALIC iterations")
		("p_weight_spatial", po::value(&opt.weight_spatial), "weight spatial")
		("p_weight_color", po::value(&opt.weight_color), "weight color")
		("p_weight_normal", po::value(&opt.weight_normal), "weight normal")
		("save_density", po::value(&p_save_density), "wether to write dasp density to file")
	;

	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	po::notify(vm);
	if(vm.count("help")) {
		std::cerr << desc << std::endl;
		return 1;
	}

	if(!p_img.empty()) {
		p_rgbd_mode = "static";
		p_rgbd_arg = p_img;
		p_stream_max = 1;
	}

	if(p_pds_mode == "rnd") {
		opt.seed_mode = dasp::SeedModes::Random;
	}
	if(p_pds_mode == "spds") {
		opt.seed_mode = dasp::SeedModes::SimplifiedPDS;
	}
	if(p_pds_mode == "delta" || p_pds_mode == "dds") {
		opt.seed_mode = dasp::SeedModes::Delta;
	}

	std::shared_ptr<RgbdStream> stream = FactorStream(p_rgbd_mode, p_rgbd_arg);

	boost::format fn_result_fmt(p_out + "%05d");

	dasp::Superpixels superpixels;
	superpixels.opt = opt;

	boost::progress_display progress(p_stream_max);

	int frame_id = 0;
	while(stream->grab()) {
		frame_id ++;
		++progress;
		if(frame_id > p_stream_max) {
			break;
		}
		Rgbd rgbd = stream->get();

		// const std::string p_img_color = p_img + "_color.png";
		// const std::string p_img_depth = p_img + "_depth.pgm";
		// if(p_verbose) std::cout << "Reading COLOR: '" << p_img_color << "'" << std::endl;
		// if(p_verbose) std::cout << "Reading DEPTH: '" << p_img_depth << "'" << std::endl;
		// slimage::Image3ub img_color = slimage::Load3ub(p_img_color);
		// slimage::Image1ui16 img_depth =  slimage::Load1ui16(p_img_depth);

		const slimage::Image3ub& img_color = rgbd.color;
		const slimage::Image1ui16& img_depth = rgbd.depth;

		if(p_verbose) std::cout << "n_given=" << opt.count << std::flush;
		dasp::ComputeSuperpixelsIncremental(superpixels, img_color, img_depth);
		if(p_verbose) std::cout << " R=" << superpixels.opt.base_radius << " n_seeds=" << superpixels.seeds.size() << " n_final=" << superpixels.cluster.size() << std::endl;

		if(!p_out.empty()) {
			const slimage::Image1i& labels = superpixels.ComputeLabels();

			std::string fn_result = (fn_result_fmt % frame_id).str();
			if(!p_img.empty()) {
				fn_result = p_out;
			}

			// image
			{
				std::string fn_img = fn_result + ".png";
				slimage::Image3ub result = dasp::plots::PlotPoints(superpixels, dasp::plots::Color);
				dasp::plots::PlotEdges(result, labels, slimage::Pixel3ub{{255,0,0}}, 1);
				if(p_verbose) std::cout << "Writing superpixel image to '" << fn_img << "'" << std::endl;
				slimage::Save(result, fn_img);
			}

			// clusters
			{
				std::string fn_clusters = fn_result + "_clusters.tsv";
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
				std::string fn_labels = fn_result + "_labels.tsv";
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

			// save density
			if(p_save_density) {
				std::string fn_density = fn_result + "_density.tsv";
				density::SaveDensity(fn_density, superpixels.density);
				if(p_verbose) std::cout << "Wrote point density to file '" << fn_density << "'." << std::endl;

				std::vector<Eigen::Vector2f> seeds(superpixels.cluster.size());
				std::transform(superpixels.cluster.begin(), superpixels.cluster.end(), seeds.begin(),
					[](const dasp::Cluster& c) {
						return Eigen::Vector2f{ c.center.px, c.center.py };
					});
				Eigen::MatrixXf approx = density::PointDensity(seeds, superpixels.density);
				std::string fn_point_density = fn_result + "_point_density.tsv";
				density::SaveDensity(fn_point_density, approx);
				if(p_verbose) std::cout << "Wrote point density to file '" << fn_point_density << "'." << std::endl;
			}
		}
	}

	return 0;
}
