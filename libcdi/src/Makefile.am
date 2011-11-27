## Process this file with automake to produce Makefile.in
#
EXTRA_DIST = make_cdilib mo_cdi.f90

if ENABLE_CDI_LIB
  lib_LTLIBRARIES = libcdi.la
  include_HEADERS = cdi.h cdi.inc
else
  noinst_LTLIBRARIES = libcdi.la
endif


libcdi_la_SOURCES = 	 \
	basetime.c     	 \
	basetime.h	 \
	binary.c	 \
	binary.h	 \
	calendar.c 	 \
	calendar.h	 \
	cdf.c            \
	cdf.h	 	 \
	cdf_int.c	 \
	cdf_int.h	 \
	cdi.h	 	 \
	cdi_error.c      \
	cdi_limits.h	 \
	cdi_util.c       \
	cdiFortran.c     \
	cfortran.h       \
	cgribex.h	 \
	cgribexlib.c  	 \
	datetime.h	 \
	dmemory.c      	 \
	dmemory.h	 \
	dtypes.h	 \
	error.c        	 \
	error.h	 	 \
	extra.h	 	 \
	extralib.c       \
	file.c         	 \
	file.h	 	 \
	gaussgrid.c      \
	gaussgrid.h	 \
	gribapi.c  	 \
	gribapi.h	 \
	grid.c           \
	grid.h	 	 \
	ieg.h	 	 \
	ieglib.c         \
	institution.c  	 \
	model.c        	 \
	service.h	 \
	servicelib.c     \
	stream_cdf.c     \
	stream_cdf.h	 \
	stream_cgribex.c \
	stream_cgribex.h \
	stream_ext.c     \
	stream_ext.h	 \
	stream_grb.c     \
	stream_grb.h     \
	stream_gribapi.c \
	stream_gribapi.h \
	stream_history.c \
	stream_ieg.c     \
	stream_ieg.h	 \
	stream_int.c     \
	stream_int.h	 \
	stream_record.c  \
	stream_srv.c     \
	stream_srv.h	 \
	stream_var.c     \
	swap.h	 	 \
	table.c        	 \
	table.h	 	 \
	tablepar.h	 \
	taxis.c          \
	taxis.h	         \
	timebase.c 	 \
	timebase.h	 \
	tsteps.c         \
	util.c         	 \
	varscan.c      	 \
	varscan.h        \
	version.c      	 \
	vlist.c 	 \
	vlist.h	         \
	vlist_att.c 	 \
	vlist_var.c 	 \
	zaxis.c		 \
        stream.c         \
        swap.c

#libcdi_la_CPPFLAGS  = @CPPFLAGS@
#libcdi_la_LIBADD    = @LDFLAGS@
#
cdilib.c:
	$(top_srcdir)/src/make_cdilib $(top_srcdir)/src
#
cdilib.o: cdilib.c
	$(COMPILE) -c $<

LOCALTARGETS  = pkgconfig/cdi.pc
if CREATE_ISOC
LOCALTARGETS += mo_cdi.o mo_cdi.$(FCMODEXT)
endif
#
mo_cdi.f90:
	ruby $(top_srcdir)/interfaces/f2003/bindGen.rb cdi.h $@
#
mo_cdi.o: mo_cdi.f90
	$(FC) $(FCFLAGS) -c $(FCFLAGS_f90) $<
#
mo_cdi.$(FCMODEXT): mo_cdi.f90
	$(FC) $(FCFLAGS) -c $(FCFLAGS_f90) $<

pkgconfig/cdi.pc: pkgconfig/cdi.pc.in ../config.status
	(cd .. ; ./config.status src/pkgconfig/cdi.pc)
#
all-local: $(LOCALTARGETS) 
#
CLEANFILES  = `ls *~`
CLEANFILES += cdilib.c
if CREATE_ISOC
CLEANFILES += $(top_builddir)/src/mo_cdi.$(FCMODEXT) $(top_builddir)/src/mo_cdi.o
endif

install-exec-local: pkgconfig/cdi.pc
	$(mkinstalldirs) "$(DESTDIR)$(libdir)/pkgconfig"
	$(install_sh_DATA) pkgconfig/cdi.pc "$(DESTDIR)$(libdir)/pkgconfig/cdi.pc"

uninstall-local:
	rm -f "$(DESTDIR)$(libdir)/pkgconfig/cdi.pc"
	rmdir "$(DESTDIR)$(libdir)/pkgconfig"

