AC_INIT([sad], [0.0.1], [gmarcais@oceangenomics.com])
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

# Boost
BOOST_REQUIRE
BOOST_STATIC
BOOST_STRING_ALGO
BOOST_IOSTREAMS

# Gurobi
AC_ARG_WITH([gurobi],
            [AC_HELP_STRING([--with-gurobi=<path>], [specify root of Gurobi installation])],
            [], [])
AS_IF([test "x$with_gurobi" != xno],
      [GUROBI=$with_gurobi]
      [GUROBI_CFLAGS="-I$GUROBI/include"]
      [GUROBI_LIB=`ls $GUROBI/lib/libgurobi@<:@0-9@:>@@<:@0-9@:>@.so | tail -n 1 | xargs basename -s .so | cut -b 4-`]
      [GUROBI_LIBS="-L$GUROBI/lib -l$GUROBI_LIB"]
      [SAVE_LIBS="$LIBS"]
      [LIBS="$LIBS $GUROBI_LIBS"]
      [AC_CHECK_LIB([$GUROBI_LIB], [GRBloadenv])]
      [LIBS="$SAVE_LIBS"])

AC_SUBST(GUROBI_CFLAGS)
AC_SUBST(GUROBI_LIBS)

# Other dependencies
PKG_CHECK_MODULES([EIGEN], [eigen3])
PKG_CHECK_MODULES([JELLYFISH], [jellyfish-2.0])
PKG_CHECK_MODULES([HTSLIB], [htslib])
PKG_CHECK_MODULES([GSL], [gsl])

# Extra RPATH stuff?
AC_LIB_LINKFLAGS_FROM_LIBS([JELLYFISH_LDFLAGS], [$JELLYFISH_LIBS])
AC_LIB_LINKFLAGS_FROM_LIBS([HTSLIB_LDFLAGS], [$HTSLIB_LIBS])
AC_LIB_LINKFLAGS_FROM_LIBS([GSL_LDFLAGS], [$GSL_LIBS])
AC_LIB_LINKFLAGS_FROM_LIBS([GUROBI_LDFLAGS], [$GUROBI_LIBS])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
