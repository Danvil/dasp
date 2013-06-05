#include <dasp/Superpixels.hpp>
#include <dasp/Plots.hpp>
#include <dasp/Segmentation.hpp>
#include <dasp/Neighbourhood.hpp>
#include <dasp/IO.hpp>
#include <rgbd/rgbd.hpp>
#include <density/PointDensity.hpp>
#include <Slimage/IO.hpp>
#include <Slimage/Slimage.hpp>
#include <Slimage/Gui.hpp>
#include <boost/program_options.hpp>
#include <boost/format.hpp>
#include <boost/progress.hpp>
#include <iostream>



int main(int argc, char** argv)
{
	constexpr slimage::Pixel3ub DASP_EDGE_COLOR{{255,0,0}};

	std::string p_img = "";
	std::string p_rgbd_mode = "";
	std::string p_rgbd_arg = "";
	std::string p_out = "";
	int p_verbose = 0;
	int p_stream_max = 100;
	std::string p_pds_mode = "spds";
	bool p_save_vis_dasp = true;
	bool p_save_cluster = true;
	bool p_save_labels = true;
	bool p_save_density = false;
	bool p_save_graph = false;

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
		("img", po::value(&p_img), "path to input image (must have X_color.png and X_depth.pgm) (overwrites rgbd_mode/_arg)")
		("rgbd_mode", po::value(&p_rgbd_mode), "rgbd stream mode (static, oni, live, ...)")
		("rgbd_arg", po::value(&p_rgbd_arg), "rgbd stream argument (depends on mode)")
		("stream_max", po::value(&p_stream_max), "maximum number of frames to process")
		("out", po::value(&p_out), "path to result files (also depends on save_* arguments)")
		("verbose", po::value(&p_verbose), "verbosity: 0 = silent, 1 = text only, 2 = basic images, 3 = all images")
		("p_radius", po::value(&opt.base_radius)->default_value(opt.base_radius), "superpixel radius [meters] (overwritten by p_count)")
		("p_count", po::value(&opt.count)->default_value(opt.count), "number of superpixels (set to 0 to use p_radius)")
		("p_pds_mode", po::value(&p_pds_mode)->default_value(p_pds_mode), "Poisson Disk sampling method (rnd, spds, dds)")
		("p_num_iterations", po::value(&opt.iterations)->default_value(opt.iterations), "number of DALIC iterations")
		("p_weight_spatial", po::value(&opt.weight_spatial)->default_value(opt.weight_spatial), "metric weight spatial")
		("p_weight_color", po::value(&opt.weight_color)->default_value(opt.weight_color), "metric weight color")
		("p_weight_normal", po::value(&opt.weight_normal)->default_value(opt.weight_normal), "metric weight normal")
		("save_vis_dasp", po::value(&p_save_vis_dasp)->default_value(p_save_vis_dasp), "enable to write dasp visualization image")
		("save_cluster", po::value(&p_save_cluster)->default_value(p_save_cluster), "enable to write dasp clusters")
		("save_labels", po::value(&p_save_labels)->default_value(p_save_labels), "enable to write dasp pixel cluster labels")
		("save_density", po::value(&p_save_density)->default_value(p_save_density), "enable to write dasp target and point density")
		("save_graph", po::value(&p_save_graph)->default_value(p_save_graph), "enable to write dasp graph to file")
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

	boost::shared_ptr<boost::progress_display> progress;
	bool show_cmd_progressbar = p_verbose == 0 && p_stream_max >= 1;
	if(show_cmd_progressbar) {
		progress.reset(new boost::progress_display(p_stream_max));
	}

	int frame_id = 0;
	while(stream->grab()) {
		// next frame
		frame_id ++;
		if(show_cmd_progressbar) {
			++(*progress);
		}
		if(frame_id > p_stream_max) {
			break;
		}

		// read next frame
		Rgbd rgbd = stream->get();
		const slimage::Image3ub& img_color = rgbd.color;
		const slimage::Image1ui16& img_depth = rgbd.depth;

		// compute superpixels
		if(p_verbose) std::cout << "n_given=" << opt.count << std::flush;
		dasp::ComputeSuperpixelsIncremental(superpixels, img_color, img_depth);
		if(p_verbose) std::cout << " R=" << superpixels.opt.base_radius << " n_seeds=" << superpixels.seeds.size() << " n_final=" << superpixels.cluster.size() << std::endl;

		// verbose / saving

		bool output_enabled = !p_out.empty();
		if(!output_enabled && (p_save_vis_dasp || p_save_density || p_save_graph || p_save_labels || p_save_cluster)) {
			std::cerr << "WARNING: Did not specify output filename! Disabled writing to files for all results." << std::endl;
		}

		slimage::Image1i labels;
		slimage::Image3ub vis_dasp;
		dasp::DaspGraph graph;

		bool needs_vis_dasp = (p_verbose >= 1) || (output_enabled && p_save_vis_dasp);
		bool needs_labels = needs_vis_dasp || (output_enabled && p_save_labels);
		bool needs_graph = (output_enabled && p_save_graph);

		// compute pixel labels
		if(needs_labels) {
			labels = superpixels.ComputeLabels();
		}

		// plot dasp superpixel image with edges
		if(needs_vis_dasp) {
			// slimage::Image3ub result = dasp::plots::PlotPoints(superpixels, dasp::plots::Color);
			vis_dasp = slimage::Image3ub(superpixels.points.width(), superpixels.points.height(), {{0,0,0}});
			dasp::plots::PlotClusters(vis_dasp, superpixels, dasp::plots::ClusterPoints, dasp::plots::Color);
			dasp::plots::PlotEdges(vis_dasp, labels, DASP_EDGE_COLOR, 2);
		}

		// compute dasp neighbourhood graph
		if(needs_graph) {
			graph = dasp::CreateDaspNeighbourhoodGraph(superpixels);
		}

		// show basic images
		if(p_verbose >= 1) {
			slimage::gui::Show("dasp", vis_dasp, 3);
		}

		if(output_enabled) {
			// filename for result files
			std::string fn_result = (fn_result_fmt % frame_id).str();
			if(!p_img.empty()) {
				fn_result = p_out;
			}

			// dasp visualization image
			if(p_save_vis_dasp) {
				std::string fn_img = fn_result + ".png";
				if(p_verbose) std::cout << "Writing superpixel image to '" << fn_img << "'" << std::endl;
				slimage::Save(vis_dasp, fn_img);
			}

			// clusters
			if(p_save_cluster) {
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
			if(p_save_labels) {
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

			// save graph
			if(p_save_graph) {
				std::string fn_vertices = fn_result + "_vertices.tsv";
				std::string fn_edges = fn_result + "_edges.tsv";
				dasp::SaveDaspGraph(graph, fn_vertices, fn_edges);
				if(p_verbose) std::cout << "Wrote dasp graph to files '" << fn_vertices << "' and '" << fn_edges << "'" << std::endl;
			}

		}
	}

	if(p_verbose >= 2) {
		slimage::gui::WaitForKeypress();
	}

	return 0;
}
