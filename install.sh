#!/bin/bash

# install the following libraries using bash in linux
# boost, eigen, gsl, jellyfish, htslib, spline

Dir=$(pwd)

mkdir -p ${Dir}/external

# install boost
echo "Install Boost library..."
cd ${Dir}/external
wget https://dl.bintray.com/boostorg/release/1.69.0/source/boost_1_69_0.tar.gz
tar -xzvf boost_1_69_0.tar.gz
cd ${Dir}/external/boost_1_69_0
mkdir Installation
./bootstrap.sh --prefix=${Dir}/external/boost_1_69_0/Installation
./b2 install


# install eigen
echo "Install Eigen library..."
cd ${Dir}/external
wget http://bitbucket.org/eigen/eigen/get/3.3.7.tar.gz
tar -xzvf 3.3.7.tar.gz
mv $(find . -maxdepth 1 -name "eigen*") eigen_3.3.7


# install gsl
echo "Install GSL library..."
cd ${Dir}/external
wget ftp://ftp.gnu.org/gnu/gsl/gsl-2.5.tar.gz
tar -xzvf gsl-2.5.tar.gz
cd ${Dir}/external/gsl-2.5/
mkdir Installation
./configure --prefix=${Dir}/external/gsl-2.5/Installation
make
make install


# install htslib
echo "Install htslib..."
cd ${Dir}/external
wget https://github.com/samtools/htslib/releases/download/1.9/htslib-1.9.tar.bz2
tar -xjvf htslib-1.9.tar.bz2
cd ${Dir}/external/htslib-1.9/
mkdir Installation
./configure --prefix=${Dir}/external/htslib-1.9/Installation
make
make install


# install jellyfish
echo "Install JELLYFISH..."
cd ${Dir}/external
wget https://github.com/gmarcais/Jellyfish/releases/download/v2.2.10/jellyfish-2.2.10.tar.gz
tar -xzvf jellyfish-2.2.10.tar.gz
cd ${Dir}/external/jellyfish-2.2.10/
mkdir Installation
./configure --prefix=${Dir}/external/jellyfish-2.2.10/Installation
make
make install


# download spline.h
echo "Download spline..."
cd ${Dir}/external
mkdir spline
cd spline
wget https://kluge.in-chemnitz.de/opensource/spline/spline.h