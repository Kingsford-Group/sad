ACLOCAL_AMFLAGS = -I m4

AM_CPPFLAGS = $(BOOST_CPPFLAGS) $(EIGEN_CFLAGS) $(JELLYFISH_CFLAGS) $(HTSLIB_CFLAGS) $(GSL_CFLAGS) $(GUROBI_CFLAGS)
AM_LDFLAGS = $(BOOST_IOSTREAMS_LDFLAGS) $(HTSLIB_LDFLAGS) $(JELLYFISH_LDFLAGS) $(GSL_LDFLAGS) $(GUROBI_LDFLAGS)
LIBS = $(BOOST_IOSTREAMS_LIBS)  $(HTSLIB_LIBS) $(GSL_LIBS) -lgurobi_g++5.2 $(GUROBI_LIBS) -lz -lm -fopenmp

bin_PROGRAMS = SAD readsalmonbias transcovdist readrsembias rsemobs

readsalmonbias_SOURCES = src/ReadSalmonBias.cpp
transcovdist_SOURCES = src/TransCovDist.cpp
readrsembias_SOURCES = src/ReadRSEMBias.cpp
rsemobs_SOURCES = src/RSEMobs.cpp
SAD_SOURCES = src/main.cpp src/DistTest.cpp src/DistTest.hpp	\
              src/Transcript.cpp src/Transcript.hpp		\
              src/LPReassign.cpp src/LPReassign.hpp src/IO.cpp	\
              src/IO.hpp
