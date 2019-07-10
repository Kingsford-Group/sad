# include folders
BOOST_INCLUDE_DIR = external/boost_1_69_0/Installation/include
EIGEN_INCLUDE_DIR = external/eigen_3.3.7
GSL_INCLUDE_DIR = external/gsl-2.5/Installation/include
GUROBI_INCLUDE_DIR = /opt/local/stow/gurobi800/linux64/include
JELLYFISH_INCLUDE_DIR = external/jellyfish-2.2.10/Installation/include/jellyfish-2.2.10
HTSLIB_INCLUDE_DIR = external/htslib-1.9/Installation/include
SPLINE_INCLUDE_DIR = external/spline

# lib folders
BOOST_LIB_DIR = external/boost_1_69_0/Installation/lib
GSL_LIB_DIR = external/gsl-2.5/Installation/lib
GUROBI_LIB_DIR = /opt/local/stow/gurobi800/linux64/lib
HTSLIB_LIB_DIR = external/htslib-1.9/Installation/lib

SHELL = /bin/bash
VER = $(shell ls ${GUROBI_LIB_DIR}/libgurobi* | grep -E "[0-9]+.so" | awk 'BEGIN{FS="[/.]"}{print substr($$(NF-1),4)}')

CC = gcc
CXX = g++
CXXFLAGS = -std=c++11 $(INCLUDES) -g
INCLUDES = -I $(BOOST_INCLUDE_DIR) -I $(EIGEN_INCLUDE_DIR) -I $(GSL_INCLUDE_DIR) -I $(GUROBI_INCLUDE_DIR) -I $(JELLYFISH_INCLUDE_DIR) -I $(HTSLIB_INCLUDE_DIR) -I $(SPLINE_INCLUDE_DIR)
LDADD = $(GSL_LIB_DIR)/libgslcblas.a $(GSL_LIB_DIR)/libgsl.a  $(HTSLIB_LIB_DIR)/libhts.a -L $(GUROBI_LIB_DIR)
LDLIBS = -lz -lm -fopenmp -lboost_iostreams -lgurobi_g++5.2 -l$(VER) -llzma -lbz2 -lcurl -lcrypto -lpthread
RPATH = $(GUROBI_LIB_DIR)

SRCS_SAD = src/main.cpp src/DistTest.cpp src/Transcript.cpp src/LPReassign.cpp src/IO.cpp
SRCS_BIAS = src/ReadSalmonBias.cpp
SRCS_COV = src/TransCovDist.cpp
SRCS_RSEMBIAS = src/ReadRSEMBias.cpp
SRCS_RSEMOBS = src/RSEMobs.cpp

all: bin/SAD bin/readsalmonbias bin/transcovdist bin/categorizesimulation

bin/SAD: $(subst .cpp,.o,$(SRCS_SAD))
	mkdir -p bin
	$(CXX) -o $@ $^ $(LDADD) $(LDLIBS) -Wl,-rpath,$(RPATH)

bin/readsalmonbias: $(subst .cpp,.o,$(SRCS_BIAS))
	mkdir -p bin
	$(CXX) -o $@ $^ $(LDADD) $(LDLIBS) -Wl,-rpath,$(RPATH)

bin/transcovdist: $(subst .cpp,.o,$(SRCS_COV))
	mkdir -p bin
	$(CXX) -o $@ $^ $(LDADD) $(LDLIBS) -Wl,-rpath,$(RPATH)

bin/readrsembias: $(subst .cpp,.o,$(SRCS_RSEMBIAS))
	mkdir -p bin
	$(CXX) -o $@ $^ $(LDADD) $(LDLIBS) -Wl,-rpath,$(RPATH)

bin/rsemobs: $(subst .cpp,.o,$(SRCS_RSEMOBS))
	mkdir -p bin
	$(CXX) -o $@ $^ $(LDADD) $(LDLIBS) -Wl,-rpath,$(RPATH)

clean:
	rm -f bin/* src/*.o
