/*
 * Evaluation.hpp
 *
 *  Created on: Mar 28, 2012
 *      Author: david
 */

#ifndef DASP_EVALUATION_HPP_
#define DASP_EVALUATION_HPP_

#include <Slimage/IO.hpp>
#include <Slimage/Slimage.hpp>
#include <Slimage/Convert.hpp>
#include <boost/format.hpp>
#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/progress.hpp>
#include <stdexcept>
#include <fstream>
#include <string>
#include <vector>

namespace dasp
{

	/** Writes a line to a CSV file */
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

	struct ProgramOptions
	{
		boost::program_options::options_description desc;

		ProgramOptions()
		: desc("Allowed options") {
			namespace po = boost::program_options;
			desc.add_options()
				("help", "produce help message")
				("db_path", po::value<std::string>(&db_path), "database base path")
				("image_path", po::value<std::string>(&image_path), "path to result images")
				("algo", po::value<std::string>(&algo), "algorithm name ('dasp' or 'slic') (required)")
				("image_ids", po::value<std::string>(), "list of image ids (e.g. \"1 3 4\" or \"3\" or \"2-7\") (required)")
				("counts", po::value<std::string>(), "list of superpixel counts (e.g. \"500 1000 2000\" or \"1000\"), default: default range")
			;
		}

		void parse(int argc, char** argv) {
			db_path = "/home/david/Documents/DataSets/dasp_rgbd_dataset";
			image_path = "/media/tmp/dasp";
			namespace po = boost::program_options;
			po::store(po::parse_command_line(argc, argv, desc), vm);
			po::notify(vm);
			if(vm.count("help")) {
				std::cerr << desc << std::endl;
				throw std::runtime_error("Help requested");
			}
			// prepare paths
			fmt_db_color = boost::format(db_path + "/images/%03d_color.png");
			fmt_db_depth = boost::format(db_path + "/images/%03d_depth.pgm");
			fmt_db_borders = boost::format(db_path + "/labels/%03d_bnds.png");
			fmt_db_labels = boost::format(db_path + "/labels/%03d_labels.pgm");
			fmt_db_labels_col = boost::format(db_path + "/labels/%03d.png");
			// get image ids
			if(vm.count("image_ids")) {
				image_ids = ParseIntegerList(vm["image_ids"].as<std::string>());
			}
			// get superpixel count
			if(vm.count("counts")) {
				counts = ParseIntegerList(vm["counts"].as<std::string>());
			}
			else {
				counts = {100,200,300,400,500,750,1000,1500,2000};
			}
		}

		bool has(const std::string& name) const {
			return vm.count(name) > 0;
		}

		template<typename T>
		T get(const std::string& name) const {
			return vm[name].as<T>();
		}

		const std::string& getAlgoName() const {
			return algo;
		}

		const std::vector<unsigned int>& getImageIds() const {
			return image_ids;
		}

		const std::vector<unsigned int>& getCounts() const {
			return counts;
		}

		slimage::Image3ub dbLoadColor(unsigned int id) const {
			return slimage::Load3ub((fmt_db_color % id).str());
		}

		slimage::Image1ui16 dbLoadDepth(unsigned int id) const {
			return slimage::Load1ui16((fmt_db_depth % id).str());
		}

		slimage::Image1ub dbLoadBorders(unsigned int id) const {
			return slimage::Pick<unsigned char>(slimage::Load((fmt_db_borders % id).str()), 0);
		}

		slimage::Image1i dbLoadLabels(unsigned int id) const {
			slimage::Image1i labels;
			slimage::conversion::Convert(slimage::Load1ui16((fmt_db_labels % id).str()), labels);
			return labels;
		}

		slimage::Image3ub dbLoadLabelsColor(unsigned int id) const {
			return slimage::Load3ub((fmt_db_labels_col % id).str());
		}

		std::string dbBasePath() const {
			return db_path;
		}

		std::string dbResultPath() const {
			return db_path + "/results";
		}

		std::string getImagePath() const {
			return image_path;
		}

	private:
		boost::program_options::variables_map vm;

		std::string db_path;

		std::string image_path;

		std::string algo;

		std::vector<unsigned int> image_ids;

		std::vector<unsigned int> counts;

		mutable boost::format fmt_db_color;
		mutable boost::format fmt_db_depth;
		mutable boost::format fmt_db_borders;
		mutable boost::format fmt_db_labels;
		mutable boost::format fmt_db_labels_col;

	private:
		/** Parses a list of integers
		 * @returns If containing a '-': list, else split by delimiters ' ,;'
		 */
		static std::vector<unsigned int> ParseIntegerList(const std::string& str)
		{
			if(str.empty()) {
				return {};
			}
			std::vector<unsigned int> img_ids;
			if(str.find('-') != std::string::npos) {
				// "a-b"
				std::vector<std::string> tokens;
				boost::split(tokens, str, boost::is_any_of("-"));
				unsigned int a = boost::lexical_cast<unsigned int>(tokens[0]);
				unsigned int b = boost::lexical_cast<unsigned int>(tokens[1]);
				for(unsigned int i=a; i<=b; i++) {
					img_ids.push_back(i);
				}
			}
			else {
				std::vector<std::string> tokens;
				boost::split(tokens, str, boost::is_any_of(" ,;"));
				for(std::string s : tokens) {
					img_ids.push_back(boost::lexical_cast<unsigned int>(s));
				}
			}
			return img_ids;
		}

	};

}

#endif
