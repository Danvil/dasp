/*
 * ObjectModel.hpp
 *
 *  Created on: Apr 20, 2012
 *      Author: david
 */

#ifndef OBJECTMODEL_HPP_
#define OBJECTMODEL_HPP_

#include "PointSet.hpp"

/** Models 3D information of an object
 * Supports three main operations:
 *  -- add superpixels
 *  -- thin out redundant superpixels
 *  -- split model using two(/several) sets of superpixels
 */
struct ObjectModel
{
	typedef std::vector<ObjectModel> ObjectModelGroup;

	void add(const Point& point) {
		points_.push_back(point);
	}

	void add(const PointSet& points) {
		points_.insert(points_.begin(), points.begin(), points.end());
	}

	void thin() {
		// TODO
	}

	ObjectModelGroup split(const std::vector<PointSet>& root_parts) const {
		// cheap: for each point in the model find nearest point in root_parts points and assign to this part
		std::vector<ObjectModel> parts(root_parts.size());
		for(const Point& p : points_) {
			float best_dist = 1e9;
			unsigned int best_index = 0;
			for(unsigned int i=0; i<root_parts.size(); i++) {
				for(const Point& q : root_parts[i]) {
					float dist = (p.position - q.position).squaredNorm();
					if(dist < best_dist) {
						best_index = i;
					}
				}
			}
			parts[best_index].points_.push_back(p);
		}
		return parts;
	}

	static ObjectModel Join(const ObjectModelGroup& parts) {
		ObjectModel u;
		for(const ObjectModel& om : parts) {
			u.add(om.points_);
		}
		return u;
	}

	PointSet points_;
};

typedef ObjectModel::ObjectModelGroup ObjectModelGroup;

#endif
