AC_INIT([sad], [0.0.2-beta], [gmarcais@oceangenomics.com])
AC_CANONICAL_HOST
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE([subdir-objects foreign parallel-tests color-tests])
AM_SILENT_RULES([yes])
AC_CONFIG_SRCDIR([src])
AC_CONFIG_HEADERS([config.h])
AC_PROG_LIBTOOL
AC_LIB_RPATH
PKG_PROG_PKG_CONFIG

# Change default compilation flags
AC_LANG(C++)
AC_PROG_CXX

# --enable-all-static
# Do not use libtool if building all static
AC_ARG_ENABLE([all-static],
              [AC_HELP_STRING([--enable-all-static], [create statically linked executable])])
STATIC_FLAGS=
AS_IF([test x$enable_all_static = xyes],
      [AC_SUBST([STATIC_FLAGS], [-all-static])]
      [PKG_CONFIG="$PKG_CONFIG --static"])

# Boost
BOOST_REQUIRE
BOOST_STATIC
BOOST_STRING_ALGO
BOOST_IOSTREAMS

# LP solvers
GUROBI_CHECK
PKG_CHECK_MODULES([CLP], [clp],
                  [HAVE_CLP=1]
                  [AC_DEFINE([HAVE_CLP], [1], [COIN-OR LP solver available])],
                  [HAVE_CLP=])

AS_IF([test x$HAVE_GUROBI = x -a x$HAVE_CLP = x],
      [AC_MSG_FAILURE("No supported LP solver was found. Please install Gurobi or Coin-OR")])

# Python3
PYTHON3_CHECK

# Other dependencies
PKG_CHECK_MODULES([EIGEN], [eigen3])
PKG_CHECK_MODULES([JELLYFISH], [jellyfish-2.0])
PKG_CHECK_MODULES([HTSLIB], [htslib])
PKG_CHECK_MODULES([GSL], [gsl])


# Extra RPATH stuff?
AS_IF([test x$enable_all_static != xyes],
      [AC_LIB_LINKFLAGS_FROM_LIBS([JELLYFISH_LDFLAGS], [$JELLYFISH_LIBS])]
      [AC_LIB_LINKFLAGS_FROM_LIBS([HTSLIB_LDFLAGS], [$HTSLIB_LIBS])]
      [AC_LIB_LINKFLAGS_FROM_LIBS([GSL_LDFLAGS], [$GSL_LIBS])]
      [AC_LIB_LINKFLAGS_FROM_LIBS([GUROBI_LDFLAGS], [$GUROBI_LIBS])]
      [AC_LIB_LINKFLAGS_FROM_LIBS([CLP_LDFLAGS], [$CLP_LIBS])])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
