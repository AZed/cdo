#  Process this file with autoconf to produce a configure script.

AC_INIT([cdi], [1.5.6], [http://code.zmaw.de/projects/cdi])

echo "configuring ${PACKAGE_NAME} ${PACKAGE_VERSION}"

CONFIG_ABORT=yes
AC_CONFIG_AUX_DIR([config])
AC_CONFIG_MACRO_DIR([m4])
AC_CANONICAL_HOST
AC_CANONICAL_BUILD

# Checks for username, hostname and system type
USERNAME=$LOGNAME
if test -z "$USERNAME" ; then USERNAME=$USER; fi;
if test -z "$USERNAME" ; then USERNAME="unknown"; fi;
AC_DEFINE_UNQUOTED([USER_NAME],["$USERNAME"], [User name])
AC_SUBST([USER_NAME],["$USERNAME"])
AS_IF([test -z "$HOST"],
      [HOST=unknown
       AS_IF([test -x /bin/hostname],[HOST=$(hostname)],
             [AS_IF([test -x /bin/uname],
                    [HOST=$(uname -n)])])])
AC_DEFINE_UNQUOTED([HOST_NAME],["$HOST"],[Host name])
AC_SUBST([HOST_NAME],["$HOST"])
AC_DEFINE_UNQUOTED([SYSTEM_TYPE],["$ac_cv_build"], [System type])
AC_SUBST([SYSTEM_TYPE],["$ac_cv_build"])

AM_INIT_AUTOMAKE
AM_CONFIG_HEADER([src/config.h])
AM_MAINTAINER_MODE([disable])

# Set up libtool.
AC_MSG_NOTICE([setting up libtool])
LT_INIT([pic-only])

# Check building environment
AC_CHECK_TOOL([CC],[gcc],[:])
AC_PROG_CC_C99
AC_PROG_FC
AS_IF([test -n "$FC" && test "X$FC" != "Xno"],
  [AC_FC_SRCEXT([f90])
   AC_LANG_PUSH([Fortran])
   AC_PROG_FPP
   AC_LANG_POP([Fortran])])
AC_PROG_F77
AC_PROG_RANLIB
AC_PROG_CXX
AC_CHECK_PROG(SED,sed,sed,false)
AC_CHECK_PROG(GREP,grep,grep,false)
AC_PROG_INSTALL
AC_CHECK_TOOL([AR],[ar],[:])
AC_CHECK_TOOL([CPP],[cpp],[:])
AC_CHECK_TOOL([LD],[ld],[:])
AC_CHECK_TOOL([NM],[nm],[:])
AC_CHECK_TOOL([AS],[as],[:])
AC_CHECK_TOOL([DLLTOOL],[dlltool],[:])
AC_CHECK_TOOL([OBJDUMP],[objdump],[:])
AC_CHECK_TOOL([STRIP],[strip],[:])
AC_CHECK_TOOL([RANLIB],[ranlib],[:])
#  ----------------------------------------------------------------------
# Check large file support on 32 bit systems
AC_SYS_LARGEFILE
#  ----------------------------------------------------------------------
# Checks for library functions.
AC_FUNC_MMAP
#  ----------------------------------------------------------------------
# Checks for library malloc.
AC_CHECK_LIB(malloc, malloc)
#  ----------------------------------------------------------------------
# Checks for structures.
AC_CHECK_MEMBERS([struct stat.st_blksize])
#  ----------------------------------------------------------------------
# Checks for header files
AC_CHECK_HEADERS(malloc.h)
#  ----------------------------------------------------------------------
# Checks for the availability of functions
AC_CHECK_FUNCS(mallinfo)
#  ----------------------------------------------------------------------
# Checks for the availability of ANSI-C99 functions
AC_CHECK_FUNCS(isnan)
#  ----------------------------------------------------------------------
#  Check for math library
AC_CHECK_LIB(m, floor)
#  ----------------------------------------------------------------------
# Add basic configure options
ACX_OPTIONS
#  ----------------------------------------------------------------------
# Compile with MPI support
AC_ARG_ENABLE(mpi,AS_HELP_STRING([--enable-mpi],[Compile with MPI compiler [default=no]]),enable_mpi=${enableval},enable_mpi=no)
#  ----------------------------------------------------------------------
#  Create the Fortran Interface via iso_c_binding module (Fortran 2003 Standard)
AC_ARG_ENABLE(iso-c-interface,
              AS_HELP_STRING([--enable-iso-c-interface],
                             [Create Fortran Interface via iso_c_bindings facility of F2003 [default=no]]).
              ,enable_isoc=${enableval},enable_isoc=no)
AM_CONDITIONAL([CREATE_ISOC],[test $enable_isoc = 'yes'])
# in case the Fortran interface uses the more modern F90-style interface,
# the FCFLAGS for users of the library have to include a switch to use the
# .mod file
CDI_F90_INTERFACE_FCFLAGS=''
AS_IF([test "x${enable_isoc}" = "xyes"],
  [AC_LANG_PUSH([Fortran])
   ACX_SL_FC_CHECK_MOD_PATH_FLAG
   AC_LANG_POP([Fortran])
   CDI_F90_INTERFACE_FCFLAGS="${FC_MOD_FLAG}${includedir}"])
AC_SUBST([CDI_F90_INTERFACE_FCFLAGS])
#  Check the module extension of the fortran compiler
ACX_SL_FC_MOD_SUFFIX(,[FCMODEXT=mod])
# -----------------------------------------------------------------------
# Check for SWIG - Generator for script-language bindings
AC_ARG_ENABLE(swig,
              [AS_HELP_STRING([--enable-swig],[use swig to create extra bindings [default=no] (EXPERIMENTAL)])],
              [AC_CHECK_PROG(SWIG,swig,swig)],
              [enable_swig=no])
AM_CONDITIONAL(ENABLE_SWIG,[test "x$SWIG" != "x"])
#  ----------------------------------------------------------------------
#  Create the Ruby Interface via swig
AC_ARG_ENABLE(ruby,
  [AS_HELP_STRING([--enable-ruby],[ruby language bindings [default=no] (EXPERIMENTAL)])],
  [AC_CHECK_PROG(RUBY,ruby,ruby)
   RUBY_INCLUDES="$($RUBY $srcdir/config/interface.rb)"
   save_CPPFLAGS=$CPPFLAGS
   CPPFLAGS="$CPPFLAGS $RUBY_INCLUDES"
   AC_CHECK_HEADER([ruby.h],,[enable_ruby=no
                              RUBY=
                             ])
   CPPFLAGS=$save_CPPFLAGS
   AS_IF([test "x$RUBY" != "x"],
         [AS_IF([test "x$SWIG" = "x"],
                [AC_MSG_NOTICE([Ruby bindings: Swig is disabled, pregenerated files will be used. Use '--enable-swig' for swig usage])])],
                [AS_IF([test `$SWIG -help 2>&1 | $GREP -c '\-ruby *- Generate'` = 0],
                       [AC_MSG_NOTICE([No (swig does not support -ruby option)])
                        enable_ruby=no
                       ])])],
  [enable_ruby=no])
AM_CONDITIONAL(ENABLE_RUBY,[test "x$enable_ruby" != "xno"])
AS_IF([test "x$enable_ruby" != "xno"],[AC_SUBST([ENABLE_RUBY],[true])],[AC_SUBST([ENABLE_RUBY],[false])])
#  ----------------------------------------------------------------------
#  Create the Python Interface via swig
AC_ARG_ENABLE(python,
  [AS_HELP_STRING([--enable-python],[python language bindings [default=no] (EXPERIMENTAL)])],
  [AC_CHECK_PROG(PYTHON,python,python)
   PYTHON_INCLUDES="-I$($PYTHON -c 'from distutils import sysconfig;print sysconfig.get_python_lib(1,1)' | $SED 's/lib/include/')"
   save_CPPFLAGS=$CPPFLAGS
   CPPFLAGS="$CPPFLAGS $PYTHON_INCLUDES"
   AC_CHECK_HEADER([Python.h],,[enable_python=no
                                PYTHON=
                               ])
   CPPFLAGS=$save_CPPFLAGS
   AS_IF([test "x$PYTHON" != "x"],
         [AS_IF([test "x$SWIG" = "x"],
                [AC_MSG_NOTICE([Python bindings: Swig is disabled, pregenerated files will be used. Use '--enable-swig' for swig usage])])],
                [AS_IF([test `$SWIG -help 2>&1 | $GREP -c '\-python *- Generate'` = 0],
                       [AC_MSG_NOTICE([No (swig does not support -python option)])
                        enable_python=no
                       ])])],
  [enable_python=no])
AM_CONDITIONAL(ENABLE_PYTHON, [test "x$enable_python" != "xno"])
AS_IF([test "x$enable_python" != "xno"],[AC_SUBST([ENABLE_PYTHON],[true])],[AC_SUBST([ENABLE_PYTHON],[false])])
#  ----------------------------------------------------------------------
AM_CONDITIONAL(CREATE_INTERFACES, [test  "x$enable_ruby" = "xyes" -o "x$enable_python" = "xyes"])
#  ----------------------------------------------------------------------
#  Check C / Fortran interface
ACX_CHECK_CFINT([$srcdir/src/cfortran.h])

AC_SUBST([CPPFLAGS])

AC_OUTPUT([Makefile src/Makefile interfaces/Makefile app/Makefile \
   tests/Makefile examples/Makefile cdi.settings src/pkgconfig/cdi.pc])

#  ----------------------------------------------------------------------
#  Show configuration
AC_MSG_NOTICE([CDI is configured with the following options:])
./config.status cdi.settings
cat cdi.settings

AS_IF([test "$ac_cv_prog_cc_c99" = "no"],
      [AC_MSG_NOTICE([

  Warning: The C compiler does not accept ANSI C99 source code!
])])
AC_MSG_NOTICE([

  Configuration completed.

    You can now say 'make' to compile the CDI package and 'make install' to install it afterwards.
])
