# Process this file with autoconf to produce a configure script.

AC_INIT(cdo, 1.4.0.1, Uwe.Schulzweida@zmaw.de)

CONFIG_ABORT=yes

echo "configuring ${PACKAGE_NAME} ${PACKAGE_VERSION}"

AC_ARG_VAR(AR,     [AR archive command]) 

AC_CONFIG_AUX_DIR(config)
AC_CANONICAL_TARGET

AM_CONFIG_HEADER(src/config.h)
AM_INIT_AUTOMAKE

# Checks for programs.
AC_PROG_CC([gcc c99])
AC_PROG_CC_C99
AC_PROG_RANLIB
AC_PROG_INSTALL
AC_CHECK_PROGS(AR, ar, [ar])

# Check large file support on 32 bit systems
AC_SYS_LARGEFILE

# Checks for library functions.
AC_FUNC_MMAP

# Checks for library malloc.
AC_CHECK_LIB(malloc, mallinfo)

# Checks for structures.
AC_CHECK_MEMBERS([struct stat.st_blksize])

# Checks for header files
AC_CHECK_HEADERS(sys/resource.h)
AC_CHECK_HEADERS(sys/times.h)

# Checks for the availability of ANSI-C99 functions
AC_CHECK_FUNCS(isnan)

AC_CHECK_FUNCS(getrlimit)

AC_CHECK_FUNCS(gethostname)

# Checks for compiler
COMPILER="$CC $CFLAGS"
AC_DEFINE_UNQUOTED(COMPILER, ["$COMPILER"], [Compiler])

# Check compiler version
case "$CC" in
  gcc*)   COMP_VERSION=`$CC --version | head -n 1`;;
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

AC_DEFINE_UNQUOTED(SYSTEM_TYPE, ["$ac_cv_target"], [System type])

#  ----------------------------------------------------------------------
#  Check for math library
#
AC_CHECK_LIB(m, floor)

#  ----------------------------------------------------------------------
#  Checks for pthread library
#
AC_ARG_WITH([threads],
            [AC_HELP_STRING([--with-threads],
                            [Use the pthread library [default=yes]])],,)

if test "$with_threads" = "no" ; then
  AC_MSG_CHECKING([for pthread library])
  AC_MSG_RESULT([suppressed])
else
  AC_CHECK_HEADERS(pthread.h)
  AC_CHECK_LIB(pthread, pthread_create)
fi

#  ----------------------------------------------------------------------
#  Checks for PROJ.4 library
#
AC_MSG_CHECKING([whether a location for the PROJ library was specified])
AC_ARG_WITH([proj],
            [AS_HELP_STRING([--with-proj=<directory>],
                            [Specify location of PROJ library. Configure will expect to find subdirs include and lib.])],
            [PROJDIR=$with_proj])
if test -z "$PROJDIR" ; then
  AC_MSG_RESULT([no])
else
  AC_MSG_RESULT([$PROJDIR])
  LDFLAGS="$LDFLAGS -L$PROJDIR/lib"
  CPPFLAGS="$CPPFLAGS -I$PROJDIR/include"
  AC_CHECK_LIB(proj, pj_init,, LDFLAGS="$LDFLAGS")
# AC_SUBST([PROJDIR], [$PROJDIR])
fi

#  ----------------------------------------------------------------------
#  Compile application with ZLIB library
#
AC_MSG_CHECKING([whether a location for the ZLIB library was specified])
AC_ARG_WITH([zlib],
            [AS_HELP_STRING([--with-zlib=<directory>],
                            [Specify location of ZLIB library. ZLIB is required for HDF5. Configure will expect to find subdirs include and lib.])],
            [ZLIBDIR=$with_zlib])
if test -z "$ZLIBDIR" ; then
  AC_MSG_RESULT([no])
else
  AC_MSG_RESULT([$ZLIBDIR])
  LDFLAGS="$LDFLAGS -L$ZLIBDIR/lib"
  CPPFLAGS="$CPPFLAGS -I$ZLIBDIR/include"
#  AC_CHECK_LIB(z, xxx, , LDFLAGS="$LDFLAGS")
# AC_SUBST([ZLIBDIR], [$ZLIBDIR])
fi

#  ----------------------------------------------------------------------
#  Compile application with SZLIB library
#
AC_MSG_CHECKING([whether a location for the SZLIB library was specified])
AC_ARG_WITH([szlib],
            [AS_HELP_STRING([--with-szlib=<directory>],
                            [Specify location of SZLIB library. You must specify its location if HDF5 was built with SZLIB. Configure will expect to find subdirs include and lib.])],
            [SZLIBDIR=$with_szlib])
if test -z "$SZLIBDIR" ; then
  AC_MSG_RESULT([no])
else
  AC_MSG_RESULT([$SZLIBDIR])
  LDFLAGS="$LDFLAGS -L$SZLIBDIR/lib"
  CPPFLAGS="$CPPFLAGS -I$SZLIBDIR/include"
  AC_CHECK_LIB(sz, SZ_BufftoBuffCompress, , LDFLAGS="$LDFLAGS")
# AC_SUBST([SZLIBDIR], [$SZLIBDIR])
fi

#  ----------------------------------------------------------------------
#  Compile application with JASPER library
#
AC_MSG_CHECKING([whether a location for the JASPER library was specified])
AC_ARG_WITH([jasper],
            [AS_HELP_STRING([--with-jasper=<directory>],
                            [Specify location of JASPER library. You must specify its location if GRIB_API was built with JASPER. Configure will expect to find subdirs include and lib.])],
            [JASPERDIR=$with_jasper])
if test -z "$JASPERDIR" ; then
  AC_MSG_RESULT([no])
else
  AC_MSG_RESULT([$JASPERDIR])
  LDFLAGS="$LDFLAGS -L$JASPERDIR/lib"
  CPPFLAGS="$CPPFLAGS -I$JASPERDIR/include"
  AC_CHECK_LIB(jasper, jas_init, , LDFLAGS="$LDFLAGS")
# AC_SUBST([JASPERDIR], [$JASPERDIR])
fi

#  ----------------------------------------------------------------------
#  Compile application with GRIB_API library
#
AC_MSG_CHECKING([whether a location for the GRIB_API library was specified])
AC_ARG_WITH([grib_api],
            [AS_HELP_STRING([--with-grib_api=<directory>],
                            [Specify location of GRIB_API library. Configure will expect to find subdirs include and lib.])],
            [GRIB_APIDIR=$with_grib_api])
if test -z "$GRIB_APIDIR" ; then
  AC_MSG_RESULT([no])
else
  AC_MSG_RESULT([$GRIB_APIDIR])
  LDFLAGS="$LDFLAGS -L$GRIB_APIDIR/lib"
  CPPFLAGS="$CPPFLAGS -I$GRIB_APIDIR/include"
  AC_CHECK_LIB(grib_api, grib_get_message, , LDFLAGS="$LDFLAGS")
# AC_SUBST([GRIB_APIDIR], [$GRIB_APIDIR])
fi

#  ----------------------------------------------------------------------
#  Compile application with HDF5 library
#
AC_MSG_CHECKING([whether a location for the HDF5 library was specified])
AC_ARG_WITH([hdf5],
            [AS_HELP_STRING([--with-hdf5=<directory>],
                            [Specify location of HDF5 library. Configure will expect to find subdirs include and lib.])],
            [HDF5DIR=$with_hdf5])
if test -z "$HDF5DIR" ; then
  AC_MSG_RESULT([no])
else
  AC_MSG_RESULT([$HDF5DIR])

  LDFLAGS="$LDFLAGS -L$HDF5DIR/lib"
  CPPFLAGS="$CPPFLAGS -I$HDF5DIR/include"
  AC_CHECK_LIB([hdf5], [H5Fopen], [], [nc_hdf5_lib_missing=yes], [-lz])
# AC_SUBST([HDF5DIR], [$HDF5DIR])
  if test "x$nc_hdf5_lib_missing" = xyes; then
    AC_MSG_ERROR([Can't find the hdf5 header and/or library, or can't link a program \
    to it. HDF5 must be built with zlib; the location of zlib must be specified for HDF5 with the \
    --with-zlib option. If HDF5 was also built with szlib, then the location of szlib must also be \
    specified with the --with-szlib option.])
  else
    LIBS="${LIBS} -lz"
# Check for the HDF5 High Level library.
    AC_CHECK_LIB([hdf5_hl], [H5DSis_scale], [], [nc_hdf5_hl_lib_missing=yes])
#    if test "x$nc_hdf5_hl_lib_missing" = xyes; then
#      AC_MSG_ERROR([Can't find the hdf5 high-level library. The high-level library is required for netCDF-4.])
#    fi
  fi

fi

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
#  Is the netcdf present? It has a header file `netcdf.h' and a library
#  `-lnetcdf' and their locations might be specified with the `--with-netcdf'
#  command-line switch. The value is an include path and/or a library path.
#  If the library path is specified then it must be preceded by a comma.
#
AC_ARG_WITH([netcdf],
            [AC_HELP_STRING([--with-netcdf=<directory>],
                            [Specify location of netCDF library. Configure will expect to find subdirs include and lib.])],,
            withval=no)

case $withval in
  yes)
    HAVE_NETCDF="yes"
    AC_CHECK_HEADERS(netcdf.h)
    AC_CHECK_LIB(netcdf, nc_open,, unset HAVE_NETCDF)

    if test -z "$HAVE_NETCDF" -a -n "$CONFIG_ABORT"; then
      AC_MSG_ERROR([Can't find netcdf library, or can't link a program to it.])
    fi
    ;;
  no)
    HAVE_NETCDF="no"
    AC_MSG_CHECKING([for netcdf library])
    AC_MSG_RESULT([suppressed])
    ;;
  *)
    HAVE_NETCDF="yes"
    case "$withval" in
      *,*)
        netcdf_inc="`echo $withval |cut -f1 -d,`"
        netcdf_lib="`echo $withval |cut -f2 -d, -s`"
        ;;
      *)
        if test -n "$withval"; then
          netcdf_inc="$withval/include"
          netcdf_lib="$withval/lib"
        fi
        ;;
    esac

#   Trying to include -I/usr/include and -L/usr/lib is redundant and
#   can mess some compilers up.
    if test "X$netcdf_inc" = "X/usr/include"; then
      netcdf_inc=""
    fi
    if test "X$netcdf_lib" = "X/usr/lib"; then
      netcdf_lib=""
    fi

#   Trying to include -I/usr/local/include and -L/usr/local/lib is redundant and
#   can mess some compilers up.
    if test "X$netcdf_inc" = "X/usr/local/include"; then
      netcdf_inc=""
    fi
    if test "X$netcdf_lib" = "X/usr/local/lib"; then
      netcdf_lib=""
    fi

    saved_CPPFLAGS="$CPPFLAGS"
    saved_LDFLAGS="$LDFLAGS"

    if test -n "$netcdf_inc"; then
      CPPFLAGS="$CPPFLAGS -I$netcdf_inc"
    fi

    AC_CHECK_HEADERS(netcdf.h,, CPPFLAGS="$saved_CPPFLAGS")

    if test -n "$netcdf_lib"; then
      LDFLAGS="$LDFLAGS -L$netcdf_lib"
    fi

    if test "${enable_dap}" = "yes"; then
      AC_CHECK_LIB(nc-dap, nc_open,, LDFLAGS="$saved_LDFLAGS"; unset HAVE_NETCDF)
if test $ac_cv_lib_nc_dap_nc_open = yes; then
cat >>confdefs.h <<_ACEOF
#define HAVE_LIBNETCDF 1
_ACEOF
fi
    else
      AC_CHECK_LIB(netcdf, nc_open,, LDFLAGS="$saved_LDFLAGS"; unset HAVE_NETCDF)
    fi

    if test -z "$HAVE_NETCDF" -a -n "$CONFIG_ABORT"; then
      AC_MSG_ERROR([Can't find netcdf library, or can't link a program to it.])
    fi
    ;;
esac

#  ----------------------------------------------------------------------

AC_SUBST([CPPFLAGS])
AC_SUBST([AR])

AC_OUTPUT(Makefile src/Makefile cdo.spec)

if test "$ac_cv_prog_cc_c99" = "no" ; then
AC_MSG_NOTICE([

  Warning: The C compiler does not accept ANSI C99 source code!
])
fi

AC_MSG_NOTICE([

  Configuration completed.

    You can now say 'make' to compile the CDO package and 'make install' to install it afterwards.
])
