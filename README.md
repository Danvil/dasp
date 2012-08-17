Depth-Adaptive Superpixels (DASP)
====

DASP is a novel oversegmentation algorithm for RGB-D images. In contrast to previous approaches it uses 3D information in addition to color. DASP can be used as a preprocessing step for various computer vision applications, for example image segmentation or point cloud segmentation.

![Superpixels and segments](https://content.wuala.com/contents/Danvil/Public/dasp/dasp.jpg)

The DASP algorithm partitions the visible surface of the 3D geometry into uniformly distributed and equally sized planar patches. This results in a classic oversegmentation of pixels into depth-adaptive superpixels which correctly reflect deformation through perspective projection.

Using for example spectral graph theory, depth-adaptive superpixels are capable to produce high-quality image and point-cloud segmentations in near realtime (~2 fps). DASP outperform state-of-the-art oversegmentation and image segmentation algorithms both in quality and runtime.

Depth-adaptive superpixel can be used in 2D as well as in 3D. The following rendering visualizes depth-adaptive superpixels in 3D space.

![3D superpoints](https://content.wuala.com/contents/Danvil/Public/dasp/dasp_3d.jpg)


Publications
----
Further technical detail can be found in the following publication:

David Weikersdorfer, David Gossow, Michael Beetz. **Depth-Adaptive Superpixels** ([pdf](https://content.wuala.com/contents/Danvil/Public/dasp/weikersdorfer2012dasp.pdf)). *21-st International Conference on Patter Recognition (ICPR), 2012*.


Installation
----
1. `git clone git://github.com/Danvil/dasp.git`
2. `cd dasp; mkdir build; cd build`
3. `cmake ..`
4. `make`

Required libraries and programs:

* Boost
Need version 1.46.1 or higher.
sudo apt-get install libboost-all-dev

* Eigen 3
sudo apt-get install libeigen3-dev

* Qt 4
sudo apt-get install libqt4-dev

* OpenCv
See http://opencv.willowgarage.com/wiki/

* OpenNI
See https://github.com/OpenNI/OpenNI
OpenNI only works for me under Ubuntu if I edit the following file
/usr/include/ni/XnPlatform.h
and replace line 74 with
#include "Linux-x86/XnPlatformLinux-x86.h"

* Misc other deps
sudo apt-get install g++ build-essentials cmake libloki-dev libloki0.1.7 libglew1.6-dev libxmu-dev
sudo apt-get install libarpack++2-dev libsuperlu3-dev


Getting started
----
Coming soon...


Dataset
----
The RGBD dataset used for this paper consists of 11 RGBD images with annotated ground truth data. All images have the resolution of 640x480 pixels and were recorded with the Microsoft Kinect sensor. Ground truth annotations were created with the [Interactive Image Segmentation Tool](http://www.eecs.berkeley.edu/Research/Projects/CS/vision/grouping/resources.html) with additional manual processing.

![001 images](https://content.wuala.com/contents/Danvil/Public/dasp/001_montage.jpg)

For each scene (ID = 001 - 011) there are the following images in the dataset:
* ID_color.png - an RGB image of the measured color values in the common PNG image format.
* ID_depth.pgm - the measured depth value stored as an uncompressed 16-bit greyscale PGM image. See the [Netpbm wikipedia article](http://en.wikipedia.org/wiki/Netpbm_format#PGM_example) for specifications.
* ID_labels.pgm - segment integer labels stored as a 16-bit greyscale PGM image.
* ID_bnds.png - segment boundaries stored as an PGM image.
* ID.png - same as ID_labels.pgm but stored as a color PNG image.

The full dataset can be downloaded [here](https://content.wuala.com/contents/Danvil/Public/dasp/dasp_rgbd_dataset.7z).

Links
----
* [David Weikersdorfer @ Technical University Munich, Germany](http://ias.cs.tum.edu/people/weikersdorfer)
* [Depth-adaptive superpixels, ICPR 2012](https://content.wuala.com/contents/Danvil/Public/dasp/weikersdorfer2012dasp.pdf)
* [DASP RGBD dataset with annotated ground truth](https://content.wuala.com/contents/Danvil/Public/dasp/dasp_rgbd_dataset.7z)
* [Interactive Image Segmentation Tool](http://www.eecs.berkeley.edu/Research/Projects/CS/vision/grouping/resources.html)