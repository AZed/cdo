# Process this file with autoconf to produce a configure script.

AC_INIT([cdo], [1.6.3], [http://code.zmaw.de/projects/cdo])

CONFIG_ABORT=yes
AC_CONFIG_AUX_DIR(config)
AC_CONFIG_MACRO_DIR([m4])
AC_CANONICAL_HOST
AC_CANONICAL_BUILD

AM_INIT_AUTOMAKE
AC_CONFIG_HEADERS([src/config.h])
AM_MAINTAINER_MODE([disable])

# Set up libtool.
AC_MSG_NOTICE([setting up libtool])
LT_INIT

# Check building environment
AC_CHECK_TOOL([CC],[gcc],[:])
AC_PROG_CC_C99
AM_PROG_CC_C_O
AC_CHECK_TOOL([AR],[ar],[:])
AC_CHECK_TOOL([CPP],[cpp],[:])
AC_CHECK_TOOL([LD],[ld],[:])
AC_CHECK_TOOL([NM],[nm],[:])
AC_CHECK_TOOL([AR],[ar],[:])
AC_CHECK_TOOL([AS],[as],[:])
AC_CHECK_TOOL([DLLTOOL],[dlltool],[:])
AC_CHECK_TOOL([OBJDUMP],[objdump],[:])
AC_CHECK_TOOL([STRIP],[strip],[:])
AC_CHECK_TOOL([RANLIB],[ranlib],[:])
AC_PROG_INSTALL
AC_PROG_CXX

FC=no
export FC
#  ----------------------------------------------------------------------
# Check large file support on 32 bit system
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
AC_CHECK_HEADERS(sys/resource.h)
AC_CHECK_HEADERS(sys/times.h)
AC_CHECK_HEADERS(malloc.h)
AC_CHECK_HEADERS(glob.h)
#  ----------------------------------------------------------------------
# Checks for the availability of functions
AC_CHECK_FUNCS(mallinfo)
#  ----------------------------------------------------------------------
# Checks for the availability of ANSI-C99 functions
AC_CHECK_DECLS([isnan],,,[AC_INCLUDES_DEFAULT
@%:@include <math.h>])
AC_CHECK_FUNCS(getrlimit)
AC_CHECK_FUNCS(gethostname)

# Check compiler version
case "$CC" in
  pgcc*)  COMP_VERSION=`$CC -V | head -2 | tail -n 1`;;
  *gcc*)  COMP_VERSION=`$CC --version | head -n 1`;;
  g++*)   COMP_VERSION=`$CC --version | head -n 1`;;
  clang*) COMP_VERSION=`$CC --version | head -n 1`;;
  icc*)   COMP_VERSION=`$CC --version | head -n 1`;;
  sxc*)   COMP_VERSION=`$CC -V 2>&1   | tail -n 1`;;
  xlc*)   COMP_VERSION=`$CC -qversion 2>&1   | head -n 1`;;
  *)      COMP_VERSION=`$CC -V 2>&1   | head -n 1` | grep -v error;;
esac

if test -z "$COMP_VERSION" ; then COMP_VERSION="unknown"; fi;
AC_DEFINE_UNQUOTED(COMP_VERSION, ["$COMP_VERSION"], [Compiler version])

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
#  ----------------------------------------------------------------------
#  Check for math library
AC_CHECK_LIB([m],[floor])
#  ----------------------------------------------------------------------
#  Enable GRIB support
AC_MSG_CHECKING([for DATA support])
AC_ARG_ENABLE([data],
              [AS_HELP_STRING([--enable-data],[DATA support [default=yes]])],
              [AS_IF([test "x$enable_data" != 'xno'],
                     [AC_DEFINE(ENABLE_DATA, [1], [Define to 1 for DATA support])
                      enable_grib=yes])],
              [AC_DEFINE(ENABLE_DATA, [1], [Define to 1 for DATA support])
               enable_data=yes])
AC_MSG_RESULT([$enable_data])
AC_SUBST([ENABLE_DATA],[$enable_data])
#  ----------------------------------------------------------------------
# Add configure options
ACX_OPTIONS
#  ----------------------------------------------------------------------
AC_SUBST([CPPFLAGS])
AC_SUBST([FCFLAGS])
AC_SUBST([CLDFLAGS])
AC_SUBST([CLIBS])
AC_SUBST([INCLUDES])
AC_SUBST([AR])

# Checks for compiler
COMPILER="$CC $CFLAGS"
AC_DEFINE_UNQUOTED(COMPILER, ["$COMPILER"], [Compiler])

AC_OUTPUT(Makefile src/Makefile contrib/Makefile test/Makefile test/data/Makefile cdo.spec cdo.settings)

#  ----------------------------------------------------------------------
#  Show configuration
AC_MSG_NOTICE([CDO is configured with the following options:])
cat cdo.settings

AS_IF([test "$ac_cv_prog_cc_c99" = "no"],
      [AC_MSG_NOTICE([

  Warning: The C compiler does not accept ANSI C99 source code!
])])
AC_MSG_NOTICE([

  Configuration completed.

    You can now say 'make' to compile the CDO package and 'make install' to install it afterwards.
])
