dnl we override this for autoconf 2.64 since some flags needed for
dnl NEC SUPER-UX are not part of the original version
dnl it's meant to stay only as long as 2.64 is still reasonably current
dnl the original Copyright notice is as follows:
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
# 02110-1301, USA.
#
# As a special exception, the Free Software Foundation gives unlimited
# permission to copy, distribute and modify the configure scripts that
# are the output of Autoconf.  You need not follow the terms of the GNU
# General Public License when using or distributing such scripts, even
# though portions of the text of Autoconf appear in them.  The GNU
# General Public License (GPL) does govern all other use of the material
# that constitutes the Autoconf program.
#
# Certain portions of the Autoconf source text are designed to be copied
# (in certain cases, depending on the input) into the output of
# Autoconf.  We call these the "data" portions.  The rest of the Autoconf
# source text consists of comments plus executable code that decides which
# of the data portions to output in any given case.  We call these
# comments and executable code the "non-data" portions.  Autoconf never
# copies any of the non-data portions into its output.
#
# This special exception to the GPL applies to versions of Autoconf
# released by the Free Software Foundation.  When you make and
# distribute a modified version of Autoconf, you may extend this special
# exception to the GPL to apply to your modified version as well, *unless*
# your modified version has the potential to copy into its output some
# of the text that was the non-data portion of the version that you started
# with.  (In other words, unless your change moves or copies text from
# the non-data portions to the data portions.)  If your modification has
# such potential, you must delete any notice of this special exception
# to the GPL from your modified version.
# _AC_FC_LIBRARY_LDFLAGS
# ----------------------
#
# Determine the linker flags (e.g. "-L" and "-l") for the Fortran
# intrinsic and runtime libraries that are required to successfully
# link a Fortran program or shared library.  The output variable
# FLIBS/FCLIBS is set to these flags.
#
# This macro is intended to be used in those situations when it is
# necessary to mix, e.g. C++ and Fortran, source code into a single
# program or shared library.
#
# For example, if object files from a C++ and Fortran compiler must
# be linked together, then the C++ compiler/linker must be used for
# linking (since special C++-ish things need to happen at link time
# like calling global constructors, instantiating templates, enabling
# exception support, etc.).
#
# However, the Fortran intrinsic and runtime libraries must be
# linked in as well, but the C++ compiler/linker doesn't know how to
# add these Fortran libraries.  Hence, the macro
# "AC_F77_LIBRARY_LDFLAGS" was created to determine these Fortran
# libraries.
#
# This macro was packaged in its current form by Matthew D. Langston.
# However, nearly all of this macro came from the "OCTAVE_FLIBS" macro
# in "octave-2.0.13/aclocal.m4", and full credit should go to John
# W. Eaton for writing this extremely useful macro.  Thank you John.
AC_DEFUN([_AC_FC_LIBRARY_LDFLAGS],
[_AC_FORTRAN_ASSERT()dnl
_AC_PROG_FC_V
AC_CACHE_CHECK([for _AC_LANG libraries of $[]_AC_FC[]], ac_cv_[]_AC_LANG_ABBREV[]_libs,
[if test "x$[]_AC_LANG_PREFIX[]LIBS" != "x"; then
  ac_cv_[]_AC_LANG_ABBREV[]_libs="$[]_AC_LANG_PREFIX[]LIBS" # Let the user override the test.
else

_AC_PROG_FC_V_OUTPUT

ac_cv_[]_AC_LANG_ABBREV[]_libs=

# Save positional arguments (if any)
ac_save_positional="$[@]"

set X $ac_[]_AC_LANG_ABBREV[]_v_output
while test $[@%:@] != 1; do
  shift
  ac_arg=$[1]
  case $ac_arg in
        [[\\/]]*.a | ?:[[\\/]]*.a)
          _AC_LIST_MEMBER_IF($ac_arg, $ac_cv_[]_AC_LANG_ABBREV[]_libs, ,
              ac_cv_[]_AC_LANG_ABBREV[]_libs="$ac_cv_[]_AC_LANG_ABBREV[]_libs $ac_arg")
          ;;
        -bI:*)
          _AC_LIST_MEMBER_IF($ac_arg, $ac_cv_[]_AC_LANG_ABBREV[]_libs, ,
             [_AC_LINKER_OPTION([$ac_arg], ac_cv_[]_AC_LANG_ABBREV[]_libs)])
          ;;
          # Ignore these flags.
        -lang* | -lcrt*.o | -lc | -lgcc* | -lSystem | -libmil | -LANG:=* | -LIST:* | -LNO:*)
          ;;
        -lkernel32)
          test x"$CYGWIN" != xyes && ac_cv_[]_AC_LANG_ABBREV[]_libs="$ac_cv_[]_AC_LANG_ABBREV[]_libs $ac_arg"
          ;;
        -[[LRuYz]])
          # These flags, when seen by themselves, take an argument.
          # We remove the space between option and argument and re-iterate
          # unless we find an empty arg or a new option (starting with -)
	  case $[2] in
	     "" | -*);;
	     *)
		ac_arg="$ac_arg$[2]"
		shift; shift
		set X $ac_arg "$[@]"
		;;
	  esac
          ;;
        -YP,*)
          for ac_j in `echo $ac_arg | sed -e 's/-YP,/-L/;s/:/ -L/g'`; do
            _AC_LIST_MEMBER_IF($ac_j, $ac_cv_[]_AC_LANG_ABBREV[]_libs, ,
                               [ac_arg="$ac_arg $ac_j"
                               ac_cv_[]_AC_LANG_ABBREV[]_libs="$ac_cv_[]_AC_LANG_ABBREV[]_libs $ac_j"])
          done
          ;;
        -[[lLR]]*)
          # -Rn for n in 0-5 controls the fortran compiler messages and
          # not the link editor on NEC SUPER-UX
          case "$host" in
            *-nec-superux)
              case $ac_arg in
                -R[[012345]])
                  ;;
                *)
                  _AC_LIST_MEMBER_IF($ac_arg, $ac_cv_[]_AC_LANG_ABBREV[]_libs, ,
                             ac_cv_[]_AC_LANG_ABBREV[]_libs="$ac_cv_[]_AC_LANG_ABBREV[]_libs $ac_arg")
              esac
              ;;
            *)
              _AC_LIST_MEMBER_IF($ac_arg, $ac_cv_[]_AC_LANG_ABBREV[]_libs, ,
                             ac_cv_[]_AC_LANG_ABBREV[]_libs="$ac_cv_[]_AC_LANG_ABBREV[]_libs $ac_arg")
              ;;
          esac
          ;;
        -[B]*)
          ac_cv_[]_AC_LANG_ABBREV[]_libs="$ac_cv_[]_AC_LANG_ABBREV[]_libs $ac_arg"
          ;;
        -ftrace | -Popenmp)
          case "$host" in
            *-nec-superux)
              _AC_LIST_MEMBER_IF($ac_arg, $ac_cv_[]_AC_LANG_ABBREV[]_libs, ,
                 ac_cv_[]_AC_LANG_ABBREV[]_libs="$ac_cv_[]_AC_LANG_ABBREV[]_libs $ac_arg")
              ;;
          esac  
          ;;
        -q32 | -q64)
          case "$host" in
            *-ibm-aix*)
              _AC_LIST_MEMBER_IF($ac_arg, $ac_cv_[]_AC_LANG_ABBREV[]_libs, ,
                 ac_cv_[]_AC_LANG_ABBREV[]_libs="$ac_cv_[]_AC_LANG_ABBREV[]_libs $ac_arg")
              ;;
          esac
          ;;
	-zallextract*| -zdefaultextract)
	  ac_cv_[]_AC_LANG_ABBREV[]_libs="$ac_cv_[]_AC_LANG_ABBREV[]_libs $ac_arg"
	  ;;
          # Ignore everything else.
  esac
done
# restore positional arguments
set X $ac_save_positional; shift

# We only consider "LD_RUN_PATH" on Solaris systems.  If this is seen,
# then we insist that the "run path" must be an absolute path (i.e. it
# must begin with a "/").
case `(uname -sr) 2>/dev/null` in
   "SunOS 5"*)
      ac_ld_run_path=`echo $ac_[]_AC_LANG_ABBREV[]_v_output |
                        sed -n 's,^.*LD_RUN_PATH *= *\(/[[^ ]]*\).*$,-R\1,p'`
      test "x$ac_ld_run_path" != x &&
        _AC_LINKER_OPTION([$ac_ld_run_path], ac_cv_[]_AC_LANG_ABBREV[]_libs)
      ;;
esac
fi # test "x$[]_AC_LANG_PREFIX[]LIBS" = "x"
])
[]_AC_LANG_PREFIX[]LIBS="$ac_cv_[]_AC_LANG_ABBREV[]_libs"
AC_SUBST([]_AC_LANG_PREFIX[]LIBS)
])# _AC_FC_LIBRARY_LDFLAGS
dnl
dnl Local Variables:
dnl mode: autoconf
dnl End:
