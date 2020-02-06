# python3.m4: Check for Python3

AC_DEFUN([PYTHON3_CHECK],
[
AC_ARG_VAR([PYTHON3], [Path to python3 executable])
AS_IF([test x$PYTHON3 = x],
      [AC_PATH_PROGS([PYTHON3], [python3 python])],
      [PYTHON3=$(realpath -s $PYTHON3)])

AC_MSG_CHECKING([Python version 3])
AS_IF([test x$PYTHON3 = x],
      [AC_MSG_RESULT([no])]
      [AC_MSG_FAILURE([No python3 executable found])])

AS_IF([$PYTHON3 --version 2>&1 | grep -q "^Python 3"],
      [AC_MSG_RESULT([yes])],
      [AC_MSG_RESULT([no])]
      [AC_MSG_FAILURE([Python executable not version 3])])

AC_SUBST([PYTHON3])
])
