noinst_PROGRAMS=collectData collectDataNStreams

AM_CPPFLAGS = -I$(top_srcdir)/src $(YAXT_CFLAGS)
AM_FCFLAGS  = $(FPP_INCOPT)$(top_srcdir)/src

if USE_MPI
noinst_PROGRAMS+=compareResourcesArray
AM_FCFLAGS += $(FPP_DEFOPT)USE_MPI
endif

if USE_FC
 noinst_PROGRAMS+=collectData2003
endif

collectData_SOURCES=collectData.c

if USE_MPI
LDADD=$(top_builddir)/src/libcdipio.la
else
LDADD=
endif
LDADD+=$(top_builddir)/src/libcdi.la

collectData2003_SOURCES=collectData2003.F90

collectDataNStreams_SOURCES=collectDataNStreams.c

compareResourcesArray_SOURCES=compareResourcesArray.c
