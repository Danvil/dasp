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

Dasp was tested under Ubuntu 11.10 and Ubuntu 12.04.

Dasp uses C++11 and requires at least GCC 4.6. Due to the poor support of the new C++ standard by Microsoft, it is probably not possible to use dasp with MSVC.

### Requirements

[Download slimage](https://content.wuala.com/contents/Danvil/Public/dasp/slimage.tar.gz) and unzip for example into the dasp root directory.

* [Boost](http://www.boost.org/) 1.46.1 or higher: `sudo apt-get install libboost-all-dev`
* [Eigen](http://eigen.tuxfamily.org) 3.x: `sudo apt-get install libeigen3-dev`
* [Qt](http://qt.nokia.com/) 4.x: `sudo apt-get install libqt4-dev`
* arpack and superlu: `sudo apt-get install libarpack++2-dev libsuperlu3-dev`
* Build essentials: `sudo apt-get install g++ build-essential cmake cmake-qt-gui`
* Misc: `sudo apt-get install libglew1.6-dev libxmu-dev`

All apt-get dependencies in one line: *sudo apt-get install libboost-all-dev libeigen3-dev libqt4-dev libarpack++2-dev libsuperlu3-dev g++ build-essential cmake cmake-qt-gui libglew1.6-dev libxmu-dev*

### Installation Instructions

1. `git clone git://github.com/Danvil/dasp.git`
2. `cd dasp; mkdir build; cd build`
3. `cmake ..`
4. `make`
5. `dasp_gui/dasp_gui` to run the Qt gui for dasp

### cmake flags

* EIGEN3_INCLUDE_DIR - Set this to the base include directory of eigen3 (e.g. '/usr/include/eigen3')
* CMAKE_BUILT_TYPE - Set to Release to compile with optimizations and get a huge speed increase
* DASP_HAS_OPENNI and OPENNI_INCLUDE_DIR - See *Kinect Live Mode* below

### Kinect Live Mode

Required if you want to process data from the Kinect in the live mode.

Download and install [OpenNI](https://github.com/OpenNI/OpenNI) and the Microsoft Kinect driver.

To enable OpenNI you have to enable the CMake flag `DASP_HAS_OPENNI` and set the CMake variable `OPENNI_INCLUDE_DIR` to the OpenNI include directory (normally `/path/to/OpenNI/Include`).


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


Issues
----

#### Compiler errors with arpack++

`/usr/include/arpack++/arrssym.h:278:7: error: ‘class ARrcSymStdEig<double>’ has no member named ‘EigVecp’`

There seems to be a bug in arpack++. Open the file `/usr/include/arpack++/arrssym.h` and go to the function `template<class ARFLOAT> int ARrcSymStdEig<ARFLOAT>::EigenValVectors(ARFLOAT* &EigVecp, ARFLOAT* &EigValp, bool ischur)`. In this function replace all occurrences of `this->EigVecp` with `EigVecp` and all occurrences of `this->EigValp` with `EigValp`.

#### Compiler errors with OpenNI

There seems to be a bug that the operating system platform is not correctly identified.

OpenNI only works for me under Ubuntu/Linux, if I edit the file `/usr/include/ni/XnPlatform.h` and replace line 74 with `#include "Linux-x86/XnPlatformLinux-x86.h"`


Links
----
* [David Weikersdorfer @ Technical University Munich, Germany](http://ias.cs.tum.edu/people/weikersdorfer)
* [Depth-adaptive superpixels, ICPR 2012](https://content.wuala.com/contents/Danvil/Public/dasp/weikersdorfer2012dasp.pdf)
* [DASP RGBD dataset with annotated ground truth](https://content.wuala.com/contents/Danvil/Public/dasp/dasp_rgbd_dataset.7z)
* [Interactive Image Segmentation Tool](http://www.eecs.berkeley.edu/Research/Projects/CS/vision/grouping/resources.html)