/*
 * register.cpp
 *
 *  Created on: Apr 10, 2012
 *      Author: david
 */

#include "DaspRegistration.hpp"
#include "KinectGrabberObject.hpp"
#include <dasp/Plots.hpp>
#include <Danvil/SimpleEngine.h>
#include <Slimage/Slimage.hpp>
#include <Slimage/IO.hpp>
#include <Slimage/Gui.hpp>
#include <Eigen/Dense>
#include <boost/bind.hpp>
#include <boost/format.hpp>
#include <boost/math/constants/constants.hpp>
#include <iostream>
#include <vector>

void RenderCorrespondence(const Pairings& pairings, const PointSet& pnts_source, const PointSet& pnts_target, const std::vector<slimage::Pixel3ub>& colors)
{
	assert(colors.size() == 0 || colors.size() == pairings.pairings_.size());
	glBegin(GL_LINES);
	for(unsigned int i=0; i<pairings.pairings_.size(); i++) {
		const std::vector<Pairings::Pair>& u = pairings.pairings_[i];
		slimage::Pixel3ub color = (colors.size() == 0) ? slimage::Pixel3ub{{255,0,0}} : colors[i];
		for(const Pairings::Pair& p : u) {
			Eigen::Vector3f a = pnts_source[p.source_id].position;
			Eigen::Vector3f n = pnts_source[p.source_id].normal;
			Eigen::Vector3f b = pnts_target[p.target_id].position;
//			glVertex3f(a[0], a[1], a[2]);
//			glVertex3f(b[0], b[1], b[2]);
			Eigen::Vector3f c = b - (b - a).dot(n) * n;
			glColor3f(0.0f, 0.0f, 0.0f);
			glVertex3f(a[0], a[1], a[2]);
			glVertex3f(c[0], c[1], c[2]);
			glColor3ub(color[0], color[1], color[2]);
			glVertex3f(c[0], c[1], c[2]);
			glVertex3f(b[0], b[1], b[2]);
		}
	}
	glEnd();
}

class IterativeClosestPointsObject
: public Danvil::SimpleEngine::IUpdateable,
  public Danvil::SimpleEngine::IRenderable
{
public:
	static constexpr float cWait = 0.0f;

	IterativeClosestPointsObject(const boost::shared_ptr<IterativeClosestPoints>& icp) {
		icp_ = icp;
		step_ = 0;
		dt_ = 0.0;
	}

	void render() {
		throw 0;
//		// render superpixels
//		glPushMatrix();
//		glMultMatrixf(icp_->T.data());
//		dasp::plots::RenderClusters(icp_->pnts_source.superpixel, dasp::plots::Color);
//		glPopMatrix();
//		dasp::plots::RenderClusters(icp_->pnts_target.superpixel, dasp::plots::Color);
//		// render point correspondences
//		RenderCorrespondence(icp_->corresponedence, icp_->T, icp_->pnts_source, Eigen::Affine3f::Identity(), icp_->pnts_target);
	}

	void update(double dt, double current_time) {
		if(icp_->isReady()) {
			return;
		}
		dt_ += dt;
		if(dt_ > cWait) {
			dt_ = 0;
			std::cout << "ICP Step " << ++step_ << std::endl;
			icp_->step();
			if(icp_->isReady()) {
				std::cout << "ICP Converged!" << std::endl;
			}
		}
	}

private:
	boost::shared_ptr<IterativeClosestPoints> icp_;
	unsigned int step_;
	double dt_;

};

class IcpBatchObject
: public Danvil::SimpleEngine::IUpdateable,
  public Danvil::SimpleEngine::IRenderable
{
public:
	IcpBatchObject(const boost::shared_ptr<DaspMultiIcp>& micp) {
		micp_ = micp;
		render_all_ = false;
	}

	bool render_all_;

	void update(double dt, double current_time) {
	}

	void onStep() {
		micp_->step();
	}

	void onNext() {
		if(micp_->isFinished()) {
			micp_->next();
		}
	}

	void render() {
//		for(int i=current_frame_, a=1; i>=0; i-=a,a++) {
//			glPushMatrix();
//			glMultMatrixf(Ts_[i].data());
//			dasp::plots::RenderClusters(dasp_points_[i].superpixel, dasp::plots::Color);
//			glPopMatrix();
//		}

//		// render superpixels
//		glPointSize(2.0f);
//		for(int i=micp_->current_frame_, a=1; i>=0; i-=a) {
//			if(!render_all_) a++;
//			PointSet u = micp_->Ts_[i].transform(micp_->dasp_points_[i].createPointSet(), micp_->partitions_[i]);
//			glBegin(GL_POINTS);
//			for(const Point& p : u.points_) {
//				glColor3fv(p.color.data());
//				glVertex3fv(p.position.data());
//			}
//			glEnd();
//		}

//		// render superpixel segments
//		{
//			glPointSize(2.0f);
//			for(int i=micp_->current_frame_, a=1; i>=0; i-=a) {
//				if(!render_all_) a++;
//				PointSet u = micp_->getFrame(i).transformPoints();
//				glBegin(GL_POINTS);
//				for(unsigned int k=0; k<u.size(); k++) {
//					const Point& p = u[k];
//					slimage::Pixel3ub color = micp_->getFrame(i).segment_colors[micp_->getFrame(i).partition.partition_uid_[k]];
//					glColor3ub(color[0], color[1], color[2]);
//					glVertex3fv(p.position.data());
//				}
//				glEnd();
//			}
//		}

		// render current superpixel segments
		{
			glPointSize(3.0f);
			PointSet u = micp_->getFrame(micp_->current_frame_).transformPoints();
			glBegin(GL_POINTS);
			for(unsigned int k=0; k<u.size(); k++) {
				const Point& p = u[k];
				slimage::Pixel3ub color = micp_->getFrame(micp_->current_frame_).segment_colors[micp_->getFrame(micp_->current_frame_).partition.partition_uid_[k]];
				glColor3ub(color[0], color[1], color[2]);
				glVertex3fv(p.position.data());
			}
			glEnd();
		}

		// render world models
		{
			glPointSize(2.0f);
			for(unsigned int k=0; k<micp_->world_.size(); k++) {
				PointSet u = micp_->world_[k].points_;
				glColor3ub(255, 255, 255);
				glBegin(GL_POINTS);
				for(const Point& p : u) {
					glVertex3fv(p.position.data());
				}
				glEnd();
			}
		}

		// renders icp pairings
		if(micp_->icp_) {
			// render point correspondences
			PointSet psrc = micp_->getFrame(micp_->current_frame_).transformPoints();
			PointSet pdst = micp_->getFrame(micp_->current_frame_-1).dasp.createPointSet();
			RenderCorrespondence(micp_->getFrame(micp_->current_frame_).pairings, psrc, pdst, micp_->getFrame(micp_->current_frame_).segment_colors);
		}

	}

private:
	boost::shared_ptr<DaspMultiIcp> micp_;
};

#include <Danvil/SimpleEngine/System/GlutSystem.h>

//#define USE_KINECT

int main(int argc, char** argv)
{
	Danvil::SimpleEngine::EnginePtr engine = Danvil::SimpleEngine::Engine::CreateDefault();

	boost::format p_fn("/home/david/Documents/DataSets/dasp_register_dataset/002_walk/kinect_%03d");

//	boost::shared_ptr<IterativeClosestPoints> icp(new IterativeClosestPoints());
//	icp->pnts_source = DaspPointSet(Rgbd((p_fn % 1).str()));
//	icp->pnts_target = DaspPointSet(Rgbd((p_fn % 3).str()));
//	engine->getScene()->addItem(icp);
//	engine->getUpdater()->addUpdateable(icp);

	std::vector<Rgbd> frames;

#ifndef USE_KINECT
	for(unsigned int i=0; i<=16; i++) {
		frames.push_back(Rgbd::Load((p_fn % i).str()));
	}
//	frames.push_back(Rgbd::Load((p_fn % 1).str()));
//	frames.push_back(Rgbd::Load((p_fn % 3).str()));
#endif

	boost::shared_ptr<DaspMultiIcp> micp(new DaspMultiIcp(frames));

	boost::shared_ptr<IcpBatchObject> icp_obj(new IcpBatchObject(micp));

#ifdef USE_KINECT
	boost::shared_ptr<KinectGrabberObject> kinect_grabber(new KinectGrabberObject());
	kinect_grabber->setNotification([micp](slimage::Image1ui16 depth, slimage::Image3ub color) {
		if(micp->isFinished()) {
			micp->addFrame(Rgbd{color, depth});
		}
	});
	icp_obj->render_all_ = true;
#endif

	engine->getScene()->addItem(icp_obj);
#ifdef USE_KINECT
	engine->getUpdater()->addUpdateable(kinect_grabber);
#endif
	engine->getUpdater()->addUpdateable(icp_obj);

	engine->getKeyboardListener()->set('s', icp_obj, &IcpBatchObject::onStep);
	engine->getKeyboardListener()->set('n', icp_obj, &IcpBatchObject::onNext);

	Danvil::SimpleEngine::GlutSystem sys(engine, "dasp register");
	return sys.main(argc, argv);
}
