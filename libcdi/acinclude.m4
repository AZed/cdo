dnl
dnl AC_CHECK_CFINT
dnl
dnl Check C / Fortran interface
dnl
AC_DEFUN([ACX_CHECK_CFINT],
  [AC_CACHE_CHECK([whether the C / Fortran interface works],[acx_cv_check_cfint],
     [AC_LANG_PUSH([C])
      AC_COMPILE_IFELSE([AC_LANG_PROGRAM([#include "$1"],[])],
        [acx_cv_check_cfint=yes],
        [acx_cv_check_cfint=no])
      AC_LANG_POP([C])])
   AS_IF([test x$acx_cv_check_cfint = xyes],
     [AC_DEFINE(HAVE_CF_INTERFACE, [1],
        [Define if C / Fortran interface cfortran.h works])])
  ])
m4_include([m4/ac_lang_program_fortran.m4])
m4_include([m4/acx_lang_fortran_check_include.m4])
m4_include([m4/acx_lang_c_check_include.m4])
dnl
dnl Local Variables:
dnl mode: autoconf
dnl End:
