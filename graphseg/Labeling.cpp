#include "Labeling.hpp"

namespace graphseg
{

void GraphLabeling::relabel()
{
	// find set of unique labels
	std::set<unsigned int> unique_labels_set(labels.begin(), labels.end());
	std::vector<unsigned int> unique_labels(unique_labels_set.begin(), unique_labels_set.end());
	num_labels = unique_labels.size();
	// create new labeling
	for(unsigned int& x : labels) {
		auto it = std::find(unique_labels.begin(), unique_labels.end(), x);
		x = it - unique_labels.begin();
	}
}

GraphLabeling GraphLabeling::CreateClean(const std::vector<unsigned int>& labels)
{
	GraphLabeling x;
	x.labels = labels;
	x.relabel();
	return x;
}

}
