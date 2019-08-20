# include folders
BOOST_INCLUDE_DIR = external/boost_1_69_0/Installation/include
EIGEN_INCLUDE_DIR = external/eigen_3.3.7
GUROBI_INCLUDE_DIR = /opt/local/stow/gurobi800/linux64/include
SPLINE_INCLUDE_DIR = external/spline

# lib folders
BOOST_LIB_DIR = external/boost_1_69_0/Installation/lib
GUROBI_LIB_DIR = /opt/local/stow/gurobi800/linux64/lib

SHELL = /bin/bash
VER = $(shell ls ${GUROBI_LIB_DIR}/libgurobi* | grep -E "[0-9]+.so" | awk 'BEGIN{FS="[/.]"}{print substr($$(NF-1),4)}')

CC = gcc
CXX = g++
CXXFLAGS = -std=c++11 $(INCLUDES) -g
INCLUDES = -I $(BOOST_INCLUDE_DIR) -I $(EIGEN_INCLUDE_DIR) $(shell pkg-config --cflags gsl) -I $(GUROBI_INCLUDE_DIR) $(shell pkg-config --cflags jellyfish-2.0) $(shell pkg-config --cflags htslib) -I $(SPLINE_INCLUDE_DIR)
LDLIBS = -lz -lm -fopenmp -L$(BOOST_LIB_DIR) -lboost_iostreams -L$(GUROBI_LIB_DIR) -lgurobi_g++5.2 -l$(VER) -lpthread $(shell pkg-config --libs gsl) $(shell pkg-config --libs htslib)
RPATH = $(BOOST_LIB_DIR):$(GUROBI_LIB_DIR):$(shell pkg-config --libs-only-L gsl | sed 's/-L//'):$(shell pkg-config --libs-only-L htslib | sed 's/-L//')

SRCS_SAD = src/main.cpp src/DistTest.cpp src/Transcript.cpp src/LPReassign.cpp src/IO.cpp
SRCS_BIAS = src/ReadSalmonBias.cpp
SRCS_COV = src/TransCovDist.cpp
SRCS_RSEMBIAS = src/ReadRSEMBias.cpp
SRCS_RSEMOBS = src/RSEMobs.cpp

all: bin/SAD bin/readsalmonbias bin/transcovdist bin/readrsembias bin/rsemobs

bin/SAD: $(subst .cpp,.o,$(SRCS_SAD))
	mkdir -p bin
	$(CXX) -o $@ $^ $(LDLIBS) -Wl,-rpath,$(RPATH)

bin/readsalmonbias: $(subst .cpp,.o,$(SRCS_BIAS))
	mkdir -p bin
	$(CXX) -o $@ $^ $(LDLIBS) -Wl,-rpath,$(RPATH)

bin/transcovdist: $(subst .cpp,.o,$(SRCS_COV))
	mkdir -p bin
	$(CXX) -o $@ $^ $(LDLIBS) -Wl,-rpath,$(RPATH)

bin/readrsembias: $(subst .cpp,.o,$(SRCS_RSEMBIAS))
	mkdir -p bin
	$(CXX) -o $@ $^ $(LDLIBS) -Wl,-rpath,$(RPATH)

bin/rsemobs: $(subst .cpp,.o,$(SRCS_RSEMOBS))
	mkdir -p bin
	$(CXX) -o $@ $^ $(LDLIBS) -Wl,-rpath,$(RPATH)

clean:
	rm -f bin/* src/*.o
