# Gurobi.m4: Locate Gurobi's library.

# Check that the Gurobi library is found. If yes, it defines
# HAVE_GUROBI. If not, noop (it is not an error).
AC_DEFUN([GUROBI_CHECK],
[AC_REQUIRE([AC_PROG_CXX])
AC_LANG_PREPROC_REQUIRE()

AC_ARG_WITH([gurobi],
            [AC_HELP_STRING([--with-gurobi=<path>], [specify root of Gurobi installation])],
            [], [])
AS_IF([test "x$with_gurobi" != x],
      [GUROBI=$with_gurobi]
      [GUROBI_CFLAGS="-I$GUROBI/include"]
      [GUROBI_LIB=`ls $GUROBI/lib/libgurobi@<:@0-9@:>@@<:@0-9@:>@.so | tail -n 1 | xargs basename -s .so | cut -b 4-`]
      [GUROBI_LIBS="-L$GUROBI/lib -l$GUROBI_LIB"])

AC_LANG_PUSH([C++])
SAVE_LIBS="$LIBS"
LIBS="$LIBS $GUROBI_LIBS"
AC_CHECK_LIB([$GUROBI_LIB], [GRBloadenv],
             [AC_DEFINE([HAVE_GUROBI], [1], [Gurobi library available])]
             [HAVE_GUROBI=1],
             [])
LIBS="$SAVE_LIBS"
AC_LANG_POP([C++])

AC_SUBST(GUROBI_CFLAGS)
AC_SUBST(GUROBI_LIBS)
])
