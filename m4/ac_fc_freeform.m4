dnl we override this for autoconf 2.64 which introduced an erroneous
dnl tab character in the definition of this test
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
#
# AC_FC_FREEFORM([ACTION-IF-SUCCESS], [ACTION-IF-FAILURE = FAILURE])
# ------------------------------------------------------------------
# Look for a compiler flag to make the Fortran (FC) compiler accept
# free-format source code, and adds it to FCFLAGS.  Call
# ACTION-IF-SUCCESS (defaults to nothing) if successful (i.e. can
# compile code using new extension) and ACTION-IF-FAILURE (defaults to
# failing with an error message) if not.  (Defined via DEFUN_ONCE to
# prevent flag from being added to FCFLAGS multiple times.)
#
# The known flags are:
#        -ffree-form: GNU g77
#                -FR: Intel compiler (icc, ecc)
#              -free: Compaq compiler (fort)
#             -qfree: IBM compiler (xlf)
# -Mfree, -Mfreeform: Portland Group compiler
#          -freeform: SGI compiler
#            -f free: Absoft Fortran
# We try to test the "more popular" flags first, by some prejudiced
# notion of popularity.
AC_DEFUN_ONCE([AC_FC_FREEFORM],
[AC_LANG_PUSH(Fortran)dnl
AC_CACHE_CHECK([for Fortran flag needed to allow free-form source],
		ac_cv_fc_freeform,
[ac_cv_fc_freeform=unknown
ac_fc_freeform_FCFLAGS_save=$FCFLAGS
for ac_flag in none -ffree-form -FR -free -qfree -Mfree -Mfreeform \
	       -freeform "-f free"
do
  test "x$ac_flag" != xnone && FCFLAGS="$ac_fc_freeform_FCFLAGS_save $ac_flag"
  AC_COMPILE_IFELSE([
  program freeform
       ! FIXME: how to best confuse non-freeform compilers?
       print *, 'Hello ', &
        'world.'
       end],
		    [ac_cv_fc_freeform=$ac_flag; break])
done
rm -f conftest.err conftest.$ac_objext conftest.$ac_ext
FCFLAGS=$ac_fc_freeform_FCFLAGS_save
])
if test "x$ac_cv_fc_freeform" = xunknown; then
  m4_default([$2],
	     [AC_MSG_ERROR([Fortran does not accept free-form source], 77)])
else
  if test "x$ac_cv_fc_freeform" != xnone; then
    FCFLAGS="$FCFLAGS $ac_cv_fc_freeform"
  fi
  $1
fi
AC_LANG_POP(Fortran)dnl
])# AC_FC_FREEFORM
dnl
dnl Local Variables:
dnl mode: autoconf
dnl license-default: ""