/*
 * Evaluation.hpp
 *
 *  Created on: Mar 28, 2012
 *      Author: david
 */

#ifndef DASP_EVALUATION_HPP_
#define DASP_EVALUATION_HPP_

#include <fstream>
#include <string>
#include <vector>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

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

/** Get image ids */
inline std::vector<unsigned int> ParseImageIds(const std::string& str)
{
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

#endif
