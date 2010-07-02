# Process this file with autoconf to produce a configure script.

AC_INIT([cdo], [1.4.5], [http://code.zmaw.de/projects/cdo])

CONFIG_ABORT=yes

AC_CONFIG_AUX_DIR(config)
AC_CANONICAL_HOST
AC_CANONICAL_BUILD

AM_INIT_AUTOMAKE
AM_CONFIG_HEADER(src/config.h)

#AC_ARG_VAR(AR,[AR archive command]) 
# Check building environment
AC_CHECK_TOOL([CC],[gcc],[:])
AC_PROG_CC_C99
#AC_PROG_RANLIB
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

# Check large file support on 32 bit system
AC_SYS_LARGEFILE

# Checks for library functions.
AC_FUNC_MMAP

# Checks for library malloc.
AC_CHECK_LIB(malloc, malloc)

# Checks for structures.
AC_CHECK_MEMBERS([struct stat.st_blksize])

# Checks for header files
AC_CHECK_HEADERS(sys/resource.h)
AC_CHECK_HEADERS(sys/times.h)
AC_CHECK_HEADERS(malloc.h)

# Checks for the availability of functions
AC_CHECK_FUNCS(mallinfo)

# Checks for the availability of ANSI-C99 functions
AC_CHECK_FUNCS(isnan)
AC_CHECK_FUNCS(getrlimit)
AC_CHECK_FUNCS(gethostname)

# Checks for compiler
COMPILER="$CC $CFLAGS"
AC_DEFINE_UNQUOTED(COMPILER, ["$COMPILER"], [Compiler])

# Check compiler version
case "$CC" in
  *gcc*)  COMP_VERSION=`$CC --version | head -n 1`;;
  g++*)   COMP_VERSION=`$CC --version | head -n 1`;;
  sxc*)   COMP_VERSION=`$CC -V 2>&1   | tail -n 1`;;
  xlc*)   COMP_VERSION=`$CC -qversion 2>&1   | head -n 1`;;
  pgcc*)  COMP_VERSION=`$CC -V | head -2 | tail -n 1`;;
  *)      COMP_VERSION=`$CC -V 2>&1   | head -n 1`;;
esac

if test -z "$COMP_VERSION" ; then COMP_VERSION="unknown"; fi;
AC_DEFINE_UNQUOTED(COMP_VERSION, ["$COMP_VERSION"], [Compiler version])

# Checks for username and hostname
USERNAME=$LOGNAME
if test -z "$USERNAME" ; then USERNAME=$USER; fi;
if test -z "$USERNAME" ; then USERNAME="unknown"; fi;

AC_DEFINE_UNQUOTED(USER_NAME, ["$USERNAME"], [User name])

if test -z "$HOST" ; then HOST="unknown"; fi;

AC_DEFINE_UNQUOTED(HOST_NAME, ["$HOST"], [Host name])

AC_DEFINE_UNQUOTED(SYSTEM_TYPE, ["$ac_cv_build"], [System type])

#  ----------------------------------------------------------------------
#  Check for math library
#
AC_CHECK_LIB(m, floor)

#  ----------------------------------------------------------------------
#  Checks for pthread library
#
AC_ARG_WITH([threads],
            [AC_HELP_STRING([--with-threads=<directory>],
                            [Use the pthread library [default=yes]])],
            [],
            [with_threads=yes])
AS_CASE([$with_threads],
        [no],[AC_MSG_CHECKING([for pthread library])
              AC_MSG_RESULT([suppressed])],
        [yes],[AC_CHECK_HEADERS(pthread.h)
               AC_CHECK_LIB(pthread, pthread_create)],
        [*],[THREADDIR=$with_threads
             LDFLAGS="$LDFLAGS -L$THREADDIR/lib"
             CPPFLAGS="$CPPFLAGS -I$THREADDIR/include"
             AC_CHECK_HEADERS(pthread.h)
             AC_CHECK_LIB(pthread, pthread_create,,LDFLAGS="$LDFLAGS")])


#  ----------------------------------------------------------------------
#  Compile application with ZLIB library
AC_ARG_WITH([zlib],
            [AS_HELP_STRING([--with-zlib=<yes|no|directory> (default=yes)],[location of ZLIB compression library (lib and include subdirs), nec. for HDF5/NETCDF4])],
            [AS_CASE(["$with_zlib"],
                     [no],[AC_MSG_CHECKING([for ZLIB library])
                           AC_MSG_RESULT([suppressed])],
                     [yes|auto],[],
                     [*],[ZLIBROOT=$with_zlib])],
            [with_zlib=yes])
AS_IF([test "x$with_zlib" != xno],
      [AC_MSG_CHECKING([for location of ZLIB library])
       HAVE_ZLIB=yes
       ACX_C_PACKAGE([zlib],[zlib.h],,,
         [HAVE_ZLIB=no
          AC_MSG_NOTICE([zlib library header not found, zlib compression not supported])],[deflate],[z],,,dnl
         [HAVE_ZLIB=no
          AC_MSG_NOTICE([cannot link zlib library, zlib compression not supported])])
       AS_IF([test $HAVE_ZLIB = yes],dnl
         [AC_DEFINE([HAVE_LIBZ],[1],[ZLIB library is present if defined to 1])
          AC_MSG_RESULT([         Using ZLIB with :  ZLIBROOT       = $ZLIBROOT])
          AC_MSG_RESULT([                            ZLIB_C_INCLUDE = $ZLIB_C_INCLUDE])
          AC_MSG_RESULT([                            ZLIB_C_LIB     = $ZLIB_C_LIB])],dnl
         [ASX_VAR_UNSET([HAVE_ZLIB])
          ASX_VAR_UNSET([ZLIB_C_INCLUDE])
          ASX_VAR_UNSET([ZLIB_C_LIB])])])
dnl ----------------------------------------------------------------------
dnl   test for SZIP compression library
dnl
AC_ARG_WITH([szlib],
            [AS_HELP_STRING([--with-szlib=<yes|no|directory> (default=yes)],[location of szlib library, optional for GRIB1 and NETCDF4 compression])],
            [AS_CASE(["$with_szlib"],
                     [no],[AC_MSG_CHECKING([for szlib library])
                           AC_MSG_RESULT([suppressed])],
                     [yes|auto],[],
                     [*],[SZLIBROOT=$with_szlib])],
            [with_szlib=no])
AS_IF([test "x$with_szlib" != xno],
      [AC_MSG_CHECKING([for location of SZLIB library])
       HAVE_SZLIB=yes
       ACX_C_PACKAGE([szlib],[szlib.h],[#ifdef HAVE_STDLIB_H
                      #include <stdlib.h>
                      #endif],,dnl
                      [HAVE_SZLIB=no
                       AC_MSG_NOTICE([szlib library header not found, szlib compression not supported])],[SZ_Decompress],[sz],,,dnl
                      [HAVE_SZLIB=no
                       AC_MSG_NOTICE([cannot link szlib library, szlib compression not supported])])
       AS_IF([test $HAVE_SZLIB = yes],dnl
             [AC_DEFINE([HAVE_LIBSZ],[1],[SZLIB library is present if defined to 1])
              AC_MSG_RESULT([         Found SZLIB with :  SZLIBROOT       = $SZLIBROOT])
              AC_MSG_RESULT([                             SZLIB_C_INCLUDE = $SZLIB_C_INCLUDE])
              AC_MSG_RESULT([                             SZLIB_C_LIB     = $SZLIB_C_LIB])],dnl
             [ASX_VAR_UNSET([HAVE_SZLIB])
              ASX_VAR_UNSET([SZLIB_C_INCLUDE])
              ASX_VAR_UNSET([SZLIB_C_LIB])])],
      [HAVE_SZLIB=yes
       ACX_C_PACKAGE([szlib],,,,,
                    [SZ_Decompress],[sz],,,dnl
                    [HAVE_SZLIB=no
                     AC_MSG_NOTICE([cannot link szlib library, szlib compression not supported])])
       AC_MSG_NOTICE([SZLIB: links to it with $SZLIB_C_LIB])])
dnl ----------------------------------------------------------------------
dnl   test for HDF5 file format library, netCDF 4.x usually requires it
dnl
AC_ARG_WITH([hdf5],
            [AS_HELP_STRING([--with-hdf5=<yes|no|directory> (default=yes)],[location of hdf5 library, nec. for import_cmsaf and NETCDF4 (hdf5 high level interface required)])],
            [AS_CASE(["$with_hdf5"],
                     [no],[AC_MSG_CHECKING([for hdf5 library])
                           AC_MSG_RESULT([suppressed])],
                     [yes|auto],[],
                     [*],[HDF5ROOT=$with_hdf5])],
            [with_hdf5=no])
AS_IF([test "x$with_hdf5" != xno],
      [AC_MSG_CHECKING([for location of HDF5 library])
       HAVE_HDF5=yes
       ACX_C_PACKAGE([hdf5],[hdf5.h],,,dnl
         [HAVE_HDF5=no
         AC_MSG_NOTICE([HDF5 library header not found, will not support netCDF 4.x])],[H5open],[hdf5],ACX_M4_GENERATE_SUBSETS([[$SZLIB_C_LIB],[$ZLIB_C_LIB]],[ ]),,dnl
         [HAVE_HDF5=no
         AC_MSG_NOTICE([cannot link HDF5 library])])
       AS_IF([test $HAVE_HDF5 = yes],dnl
         [AC_DEFINE([HAVE_LIBHDF5],[1],[HDF 5 library is present if defined to 1])
          AC_MSG_RESULT([         Using HDF5 with :  HDF5ROOT       = $HDF5ROOT])
          AC_MSG_RESULT([                            HDF5_C_INCLUDE = $HDF5_C_INCLUDE])
          AC_MSG_RESULT([                            HDF5_C_LIB     = $HDF5_C_LIB])],dnl
         [ASX_VAR_UNSET([HAVE_HDF5])
          ASX_VAR_UNSET([HDF5_C_INCLUDE])
          ASX_VAR_UNSET([HDF5_C_LIB])])],
      [HAVE_HDF5=yes
       ACX_C_PACKAGE([hdf5],,,,dnl
         [HAVE_HDF5=no
         AC_MSG_NOTICE([HDF5 library header not found, will not support netCDF 4.x])],[H5open],[hdf5],ACX_M4_GENERATE_SUBSETS([[$SZLIB_C_LIB],[$ZLIB_C_LIB]],[ ]),,dnl
         [HAVE_HDF5=no
         AC_MSG_NOTICE([cannot link HDF5 library])])])
dnl ----------------------------------------------------------------------
dnl   test for netCDF file format library
dnl
AC_ARG_WITH([netcdf],
            [AS_HELP_STRING([--with-netcdf=<yes|no|directory> (default=yes)],[location of netcdf library (lib and include subdirs)])],
            [AS_CASE(["$with_netcdf"],
                     [no],[AC_MSG_CHECKING([for netcdf library])
                           AC_MSG_RESULT([suppressed])],
                     [yes|auto],[],
                     [*],[NETCDFROOT=$with_netcdf])],
            [with_netcdf=no])
AS_IF([test "x$with_netcdf" != xno],
      [AC_MSG_CHECKING([for location of NETCDF Library])
       HAVE_NETCDF=yes
       ACX_C_PACKAGE([netcdf],[netcdf.h],,,dnl
         [HAVE_NETCDF=no
          AC_MSG_FAILURE([Required include netcdf.h not found or not compilable.])],[nc_open],[netcdf],ACX_M4_GENERATE_SUBSETS([[-lnetcdf],[-lhdf5_hl],[$HDF5_C_LIB]],[ ]),,dnl
         [HAVE_NETCDF=no
          AC_MSG_FAILURE([Cannot link with NetCDF library.])],[])
       AS_IF([test $HAVE_NETCDF = yes],dnl
         [AC_DEFINE([HAVE_LIBNETCDF],[1],[NETCDF library is present if defined to 1])
          AC_MSG_RESULT([         Using NETCDF with : NETCDFROOT       = $NETCDFROOT])
          AC_MSG_RESULT([                             NETCDF_C_INCLUDE = $NETCDF_C_INCLUDE])
          AC_MSG_RESULT([                             NETCDF_C_LIB     = $NETCDF_C_LIB])],dnl
         [ASX_VAR_UNSET([HAVE_NETCDF])
          ASX_VAR_UNSET([NETCDF_C_INCLUDE])
          ASX_VAR_UNSET([NETCDF_C_LIB])])])
#  ----------------------------------------------------------------------
#  Checks for PROJ.4 library
#
AC_ARG_WITH([proj],
            [AS_HELP_STRING([--with-proj=<yes|no|directory>],
                            [library for cartographic projections; if a directory is given, it will be used as a value for --with-proj-root])],
            [AS_CASE(["$with_proj"],
                     [no],[AC_MSG_CHECKING([for proj library])
                           AC_MSG_RESULT([suppressed])],
                     [yes|auto],[],
                     [*],[PROJROOT=$with_proj])],
            [with_proj=yes])
AS_IF([test "x$with_proj" != xno],
      [AC_MSG_CHECKING([for location of PROJ library])
       HAVE_PROJ=yes
       ACX_C_PACKAGE([proj],[projects.h],,,
                     [HAVE_PROJ=no
                      AC_MSG_NOTICE([proj library header not found, remapping not supported])],[pj_init],[proj],,,dnl
                     [HAVE_PROJ=no
                      AC_MSG_NOTICE([cannot link proj library, remapping not supported])])
       AS_IF([test $HAVE_PROJ = yes],
             [AC_DEFINE([HAVE_LIBPROJ],[1],[PROJ library is present if defined to 1])
              AC_MSG_RESULT([         Using PROJ with : PROJROOT       = $PROJROOT])
              AC_MSG_RESULT([                           PROJ_C_INCLUDE = $PROJ_C_INCLUDE])
              AC_MSG_RESULT([                           PROJ_C_LIB     = $PROJ_C_LIB])],dnl
             [ASX_VAR_UNSET([HAVE_PROJ])
              ASX_VAR_UNSET([PROJ_C_INCLUDE])
              ASX_VAR_UNSET([PROJ_C_LIB])
              AC_MSG_RESULT([not found])],
              )])
#  ----------------------------------------------------------------------
#  Compile application with JASPER library
AC_ARG_WITH([jasper],
            [AS_HELP_STRING([--with-jasper=<yes|no|directory>],
                            [library for grib2 compression; if a directory is given, it will be used as a value for --with-jasper-root])],
            [AS_CASE(["$with_jasper"],
                     [no],[AC_MSG_CHECKING([for jasper library])
                           AC_MSG_RESULT([suppressed])],
                     [yes|auto],[],
                     [*],[JASPERROOT=$with_jasper])],
            [with_jasper=yes])
AS_IF([test "x$with_jasper" != xno],
      [AC_MSG_CHECKING([for location of JASPER library])
       HAVE_JASPER=yes
       ACX_C_PACKAGE([jasper],,,,,[jas_init],[jasper],,,dnl
                     [HAVE_JASPER=no
                      AC_MSG_NOTICE([cannot link jasper library, jpeg compression for grib2 not supported])])]
       AS_IF([test $HAVE_JASPER = no],
             [ASX_VAR_UNSET([HAVE_JASPER])
              ASX_VAR_UNSET([JASPER_C_INCLUDE])
              ASX_VAR_UNSET([JASPER_C_LIB])
              AC_MSG_RESULT([not found])]))
#  ----------------------------------------------------------------------
#  Compile application with GRIB_API library
#
AC_ARG_WITH([grib-api],
            [AS_HELP_STRING([--with-grib-api=<yes|no|directory>],
                            [Specify location of GRIB_API library (GRIB2 support). Configure will expect to find subdirs include and lib.])],
            [AS_CASE(["$with_grib2"],
                     [no],[AC_MSG_CHECKING([for grib_api library])
                           AC_MSG_RESULT([suppressed])],
                     [yes|auto],[],
                     [*],[GRIB_APIROOT=$with_grib_api])],
            [with_grib_api=no])
AS_IF([test "x$with_grib_api" != xno],
      [AC_MSG_CHECKING([for location of GRIB_API library])
       HAVE_GRIB_API=yes
       ACX_C_PACKAGE([grib_api],[grib_api.h],,,
                     [HAVE_GRIB_API=no
                      AC_MSG_NOTICE([grib_api library header not found in "$GRIP_API_C_INCLUDE", remapping not supported])],[grib_get_message],[grib_api],ACX_M4_GENERATE_SUBSETS([[$JASPER_C_LIB]],[ ]),,dnl
                     [HAVE_GRIB_API=no
                      AC_MSG_NOTICE([cannot link grib_api library, remapping not supported])])
       AS_IF([test $HAVE_GRIB_API = yes],
             [AC_DEFINE([HAVE_LIBGRIB_API],[1],[GRIB_API library is present if defined to 1])
              AC_MSG_RESULT([         Using GRIB_API with : GRIB_APIROOT       = $GRIB_APIROOT])
              AC_MSG_RESULT([                               GRIB_API_C_INCLUDE = $GRIB_API_C_INCLUDE])
              AC_MSG_RESULT([                               GRIB_API_C_LIB     = $GRIB_API_C_LIB])],dnl
             [ASX_VAR_UNSET([HAVE_GRIB_API])
              ASX_VAR_UNSET([GRIB_API_C_INCLUDE])
              ASX_VAR_UNSET([GRIB_API_C_LIB])
              AC_MSG_RESULT([not found])],
              )])
#  ----------------------------------------------------------------------
#  Enable GRIB support
AC_ARG_ENABLE(grib,AS_HELP_STRING([--enable-grib],[GRIB support [[default=yes]]]),enable_grib=${enableval},enable_grib=yes)

if test "${enable_grib}" = "yes"; then
  AC_DEFINE(HAVE_LIBGRIB, [1], [Define to 1 for GRIB support])
fi

AC_MSG_CHECKING([for GRIB support])
AC_MSG_RESULT($enable_grib)

#  ----------------------------------------------------------------------
#  Compile interface with internal CGRIBEX library
#
AC_ARG_WITH([cgribex],
            [AC_HELP_STRING([--with-cgribex], [Use the GRIB library [default=yes]])],,
             with_cgribex="yes")

if test "${with_cgribex}" = "yes" ; then
  AC_DEFINE(HAVE_LIBCGRIBEX, [1], [Define to 1 for GRIB1 decoding/encoding with cgribex])
fi

AC_MSG_CHECKING([for CGRIBEX support])
AC_MSG_RESULT($with_cgribex)

#  ----------------------------------------------------------------------
#  Compile interface with internal SERVICE library
#
AC_ARG_WITH([service],
            [AC_HELP_STRING([--with-service],
                            [Use the internal SERVICE library [default=yes]])],,
             with_service="yes")

if test "${with_service}" = "yes" ; then
  AC_DEFINE(HAVE_LIBSERVICE, [1], [Define to 1 for SERVICE interface])
fi

AC_MSG_CHECKING([for SERVICE support])
AC_MSG_RESULT($with_service)

#  ----------------------------------------------------------------------
#  Compile interface with internal EXTRA library
#
AC_ARG_WITH([extra],
            [AC_HELP_STRING([--with-extra],
                            [Use the internal EXTRA library [default=yes]])],,
             with_extra="yes")

if test "${with_extra}" = "yes" ; then
  AC_DEFINE(HAVE_LIBEXTRA, [1], [Define to 1 for EXTRA interface])
fi

AC_MSG_CHECKING([for EXTRA support])
AC_MSG_RESULT($with_extra)

#  ----------------------------------------------------------------------
#  Compile interface with internal IEG library
#
AC_ARG_WITH([ieg],
            [AC_HELP_STRING([--with-ieg],
                            [Use the internal IEG library [default=yes]])],,
             with_ieg="yes")

if test "${with_ieg}" = "yes" ; then
  AC_DEFINE(HAVE_LIBIEG, [1], [Define to 1 for IEG interface])
fi

AC_MSG_CHECKING([for IEG support])
AC_MSG_RESULT($with_ieg)

#  ----------------------------------------------------------------------
# Begin DAP (aka DODS) support
# Store pre-DAP LIBS value to restore from later
LIBS_no_DAP="${LIBS}"

AC_ARG_ENABLE(dap,AS_HELP_STRING([--enable-dap],[Build OPeNDAP/DODS-enabled clients if possible [[default=no]]]),enable_dap=${enableval},enable_dap=no)

if test "${enable_dap}" = "yes"; then
   if which ncdap-config > /dev/null; then
      CPPFLAGS="${CPPFLAGS} `ncdap-config --cflags`"
      nc_dap_ver=`ncdap-config --version`
      LIBS="`ncdap-config --libs` ${LIBS}"	
   else
      AC_MSG_ERROR(couldn't find ncdap-config)     
   fi
fi

#  ----------------------------------------------------------------------
AC_SUBST([CPPFLAGS])
AC_SUBST([FCFLAGS])
AC_SUBST([CLDFLAGS])
AC_SUBST([CLIBS])
AC_SUBST([INCLUDES])
AC_SUBST([AR])

AC_OUTPUT(Makefile src/Makefile contrib/Makefile cdo.spec)

if test "$ac_cv_prog_cc_c99" = "no" ; then
AC_MSG_NOTICE([

  Warning: The C compiler does not accept ANSI C99 source code!
])
fi

AC_MSG_NOTICE([

  Configuration completed.

    You can now say 'make' to compile the CDO package and 'make install' to install it afterwards.
])
