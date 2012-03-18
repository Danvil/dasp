#include <boost/lexical_cast.hpp>
#include <SuperPoints/Superpixels.hpp>
#include <SuperPoints/Plots.hpp>
#include <Slimage/IO.hpp>
#include <Slimage/Slimage.hpp>
#include <boost/format.hpp>
#include <boost/program_options.hpp>
#include <iostream>

dasp::Clustering ComputeClustering(slimage::Image3ub color, slimage::Image1ui16 depth, unsigned int cnt)
{
	dasp::Clustering clustering;
	clustering.opt.camera = dasp::Camera{320.0f, 240.0f, 540.0f, 0.001f};
	clustering.opt.coverage = 1.7f;
	clustering.opt.weight_spatial = 1.0f;
	clustering.opt.weight_color = 2.0f;
	clustering.opt.weight_normal = 3.0f;
	clustering.opt.iterations = 50;
	clustering.opt.seed_mode = dasp::SeedModes::DepthMipmap;
	clustering.opt.base_radius = 0.012f;
	clustering.opt.count = cnt;

	clustering.CreatePoints(color, depth);

	std::vector<dasp::Seed> seeds = clustering.FindSeeds();

	clustering.ComputeSuperpixels(seeds);
	clustering.ComputeExt();

	return clustering;
}

slimage::Image3ub ComputeBoundaryImage(const dasp::Clustering& clustering)
{
	slimage::Image3ub img_bnds(clustering.width(), clustering.height());
	img_bnds.fill({{0,0,0}});
	dasp::plots::PlotEdges(img_bnds, clustering.ComputePixelLabels(), slimage::Pixel3ub{{255,255,255}}, 1);
	return img_bnds;
}

template<typename K>
void WriteCsvLine(std::ofstream& ofs, const std::vector<K>& vals) {
	if(vals.size() > 0) {
		ofs << vals[0];
		for(unsigned int i=1; i<vals.size(); i++) {
			ofs << "," << vals[i];
		}
	}
	ofs << std::endl;
}

template<typename K>
void WriteCsvLine(std::ofstream& ofs, unsigned int count, const std::vector<K>& vals) {
	ofs << count;
	for(K v : vals) {
		ofs << "," << v;
	}
	ofs << std::endl;
}

int main(int argc, char** argv)
{
	std::string p_fpath = "/home/david/Documents/DataSets/dasp_rgbd_dataset";
	unsigned int p_db_img_id = 2;
	std::string p_db_algo_name = "turbo";

	namespace po = boost::program_options;
	po::options_description desc("Allowed options");
	desc.add_options()
	    ("help", "produce help message")
	    ("db_base_path", po::value<std::string>(&p_fpath), "database mode: base path")
	    ("db_img_id", po::value<unsigned int>(&p_db_img_id), "database mode: image id")
	    ("db_algo", po::value<std::string>(&p_db_algo_name), "database mode: algorithm name")
	;
	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	po::notify(vm);
	if (vm.count("help")) {
		std::cout << desc << std::endl;
	    return 1;
	}

	std::vector<unsigned int> counts{100,200,300,400,500,750,1000,1500,2000};

	std::string id_str = (boost::format("%03d") % p_db_img_id).str();
	std::string fn_color = p_fpath + "/images/" + id_str + "_color.png";
	std::string fn_depth = p_fpath + "/images/" + id_str + "_depth.pgm";
	std::string fn_hist_thick = p_fpath + "/results/" + id_str + "_" + p_db_algo_name + "_thickness.csv";
	std::string fn_hist_discr = p_fpath + "/results/" + id_str + "_" + p_db_algo_name + "_discR.csv";

	unsigned int p_hist_cnt = 128;
	float p_max_thick = 0.02f;
	slimage::Image3ub img_col = slimage::Load3ub(fn_color);
	slimage::Image1ui16 img_dep = slimage::Load1ui16(fn_depth);
	std::ofstream of_thick(fn_hist_thick);
	std::ofstream of_disc(fn_hist_discr);

	std::vector<float> thick_header(p_hist_cnt);
	for(unsigned int i=0; i<p_hist_cnt; i++) {
		thick_header[i] = p_max_thick * static_cast<float>(i) / static_cast<float>(p_hist_cnt - 1);
	}
	WriteCsvLine(of_thick, 0, thick_header);
	std::vector<float> discr_header(p_hist_cnt);
	for(unsigned int i=0; i<p_hist_cnt; i++) {
		discr_header[i] = 0.5f + 1.5f * static_cast<float>(i) / static_cast<float>(p_hist_cnt - 1);
	}
	WriteCsvLine(of_disc, 0, discr_header);

	for(unsigned int count : counts) {
		std::cout << "N=" << count << " ..." << std::endl;
		std::string count_str = (boost::format("%d") % count).str();
		std::string fn_bnds = p_fpath + "/results/" + id_str + "_" + p_db_algo_name + "_" + count_str + ".png";
		dasp::Clustering clustering = ComputeClustering(img_col, img_dep, count);
		// boundary image
		slimage::Image3ub bnd = ComputeBoundaryImage(clustering);
		slimage::Save(bnd, fn_bnds);
//		// histograms
//		dasp::ClusterGroupInfo info = clustering.ComputeClusterGroupInfo(p_hist_cnt, p_max_thick);
//		WriteCsvLine(of_thick, count, info.hist_thickness.bins());
//		WriteCsvLine(of_disc, count, info.hist_coverage.bins());
		// write eigenvalues
		std::string fn_eigenvalues = p_fpath + "/results/" + id_str + "_" + p_db_algo_name + "_" + count_str + "_clusters.csv";
		std::ofstream of_eigenvalues(fn_eigenvalues);
		clustering.ForClustersNoReturn([&of_eigenvalues](const dasp::Cluster& c) {
			WriteCsvLine(of_eigenvalues, std::vector<float>{c.ew(0), c.ew(1), c.ew(2), c.thickness, c.circularity, c.coverage_error});
		});
	}

	return 0;
}
