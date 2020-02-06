ACLOCAL_AMFLAGS = -I m4

AM_CPPFLAGS = -std=c++11 $(BOOST_CPPFLAGS) $(EIGEN_CFLAGS) $(JELLYFISH_CFLAGS) $(HTSLIB_CFLAGS) $(GSL_CFLAGS) $(GUROBI_CFLAGS) $(CLP_CFLAGS) -fopenmp
AM_LDFLAGS = $(BOOST_IOSTREAMS_LDFLAGS) $(HTSLIB_LDFLAGS) $(JELLYFISH_LDFLAGS) $(GSL_LDFLAGS) $(GUROBI_LDFLAGS) $(CLP_LDFLAGS) $(STATIC_FLAGS) -fopenmp
LIBS = $(BOOST_IOSTREAMS_LIBS)  $(HTSLIB_LIBS) $(GSL_LIBS) $(GUROBI_LIBS) $(CLP_LIBS) -lz -lm -fopenmp

EXTRA_DIST = src/spline.h src/logtime.hpp src/massert.hpp install-dependencies.sh
CLEANFILES =

################
# C++ programs #
################
bin_PROGRAMS = SAD readsalmonbias transcovdist readrsembias rsemobs

readsalmonbias_SOURCES = src/ReadSalmonBias.cpp
transcovdist_SOURCES = src/TransCovDist.cpp
readrsembias_SOURCES = src/ReadRSEMBias.cpp
rsemobs_SOURCES = src/RSEMobs.cpp
SAD_SOURCES = src/main.cpp src/DistTest.cpp src/DistTest.hpp	\
              src/Transcript.cpp src/Transcript.hpp		\
              src/LPReassign.cpp src/LPReassign.hpp src/IO.cpp	\
              src/IO.hpp

###########
# Scripts #
###########
# For each script, generate 2 versions: one with local path, to use
# directly; and one (in .libs/) to use after installation. Similar to
# what LIBTOOLS does.

python_scripts = AdjustPValue.py plot_coverage.py		\
                 PlotDisttributionTrans.py Pred_Trans2Gene.py	\
                 SADpipe.py TranscriptClass.py
all_scripts = $(python_scripts)
libs_scripts = $(patsubst %,.libs/%,$(all_scripts))
bin_SCRIPTS = $(libs_scripts)
noinst_SCRIPTS = $(all_scripts)
EXTRA_DIST += $(patsubst %,scripts/%,$(all_scripts))
CLEANFILES += $(libs_scripts) $(all_scripts)


# Rules to generate the scripts, in .libs and in build directory.
script_inst_subst  = sed -e 's,[@]PYTHON3[@],$(PYTHON3),g' \
                       -e 's,[@]BIN_DIR[@],$(bindir),g'
script_local_subst = sed -e 's,[@]PYTHON3[@],$(PYTHON3),g' \
                       -e 's,[@]BIN_DIR[@],$(abs_builddir),g'

.libs/%.py: scripts/%.py
	@mkdir -p $(dir $@) && $(script_inst_subst) < $< > $@
%.py: scripts/%.py
	$(AM_V_GEN)$(script_local_subst) < $< > $@ && chmod a+x $@

# Extra flag for the development tree
-include $(top_srcdir)/dev.mk
