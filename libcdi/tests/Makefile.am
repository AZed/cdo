## Process this file with automake to produce Makefile.in
#
TESTS = test_grib.sh
check_PROGRAMS = test_grib
#
test_grib_SOURCES = test_grib.c
#
LDADD = -L../src -lcdi -lm
AM_LDFLAGS = -L$(top_srcdir)/src
INCLUDES   = -I$(top_srcdir)/src
#
EXTRA_DIST = $(TESTS)
#
CLEANFILES  = `ls *~ *.grb *.nc *.srv *.ext`
