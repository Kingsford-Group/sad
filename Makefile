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

CC = gcc
CXX = g++
CXXFLAGS = -std=c++11 $(INCLUDES) -g
INCLUDES = -I $(BOOST_INCLUDE_DIR) -I $(EIGEN_INCLUDE_DIR) -I $(GSL_INCLUDE_DIR) -I $(GUROBI_INCLUDE_DIR) -I $(JELLYFISH_INCLUDE_DIR) -I $(HTSLIB_INCLUDE_DIR) -I $(SPLINE_INCLUDE_DIR)
LDADD = $(GSL_LIB_DIR)/libgslcblas.a $(GSL_LIB_DIR)/libgsl.a  $(HTSLIB_LIB_DIR)/libhts.a -L $(GUROBI_LIB_DIR)
LDLIBS = -lz -lm -fopenmp -lboost_iostreams -lgurobi_c++ -lgurobi80 -llzma -lbz2 -lcurl -lcrypto -lpthread
RPATH = $(GUROBI_LIB_DIR)

SRCS_SAD = src/main.cpp src/DistTest.cpp src/Transcript.cpp src/LPReassign.cpp src/IO.cpp
SRCS_BIAS = src/ReadSalmonBias.cpp
SRCS_COV = src/TransCovDist.cpp

all: bin/SAD bin/readsalmonbias bin/transcovdist

bin/SAD: $(subst .cpp,.o,$(SRCS_SAD))
	mkdir -p bin
	$(CXX) -o $@ $^ $(LDADD) $(LDLIBS) -Wl,-rpath,$(RPATH)

bin/readsalmonbias: $(subst .cpp,.o,$(SRCS_BIAS))
	mkdir -p bin
	$(CXX) -o $@ $^ $(LDADD) $(LDLIBS) -Wl,-rpath,$(RPATH)

bin/transcovdist: $(subst .cpp,.o,$(SRCS_COV))
	mkdir -p bin
	$(CXX) -o $@ $^ $(LDADD) $(LDLIBS) -Wl,-rpath,$(RPATH)

clean:
	rm -f bin/* src/*.o
