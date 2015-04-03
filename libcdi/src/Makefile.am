## Process this file with automake to produce Makefile.in
#
EXTRA_DIST = cdilib.c mo_cdi.f90

if ENABLE_CDI_LIB
  lib_LTLIBRARIES = libcdi.la
  include_HEADERS = cdi.h cdi.inc
else
  noinst_LTLIBRARIES = libcdi.la
endif

AM_CFLAGS = $(PPM_CORE_CFLAGS) $(YAXT_CFLAGS)

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
	cdiFortran.h     \
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
	institution.h  	 \
	model.c        	 \
	model.h        	 \
	namespace.c      \
	namespace.h      \
	serialize.h	\
	serialize.c	\
	pio.c            \
	pio.h            \
	pio_comm.c       \
	pio_comm.h       \
	pio_dbuffer.c    \
	pio_impl.h 	 \
	pio_interface.c  \
	pio_interface.h	 \
	pio_mpinonb.c    \
	pio_record_send.c\
	pio_posixasynch.c\
	pio_posixfpguardsendrecv.c \
	pio_posixnonb.c  \
	pio_list_set.c	 \
	resource_unpack.h \
	resource_unpack.c \
	pio_rpc.h        \
	pio_client.c	\
	pio_client.h	\
	pio_server.c     \
	pio_server.h     \
	pio_util.c       \
	pio_util.h       \
	pio_serialize.h	\
	pio_serialize.c	\
	pio_cdf_int.c	\
	pio_cdf_int.h	\
	resource_handle.c\
	resource_handle.h\
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
	cdi_int.c        \
	cdi_int.h	 \
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
	vlist_att.h 	 \
	vlist_var.c 	 \
	vlist_var.h	 \
	zaxis.c		 \
	zaxis.h		 \
        stream.c         \
        swap.c

#libcdi_la_CPPFLAGS  = @CPPFLAGS@
libcdi_la_LIBADD    = $(PPM_CORE_LIBS) $(YAXT_LIBS)
#
#cdilib.c:
#	$(top_srcdir)/src/make_cdilib $(top_srcdir)/src
#
#cdilib.o: cdilib.c
#	$(COMPILE) -c $<

LOCALTARGETS  = 
if ENABLE_CDI_LIB
LOCALTARGETS += pkgconfig/cdi.pc
endif

if CREATE_ISOC
LOCALTARGETS += mo_cdi.o mo_cdi.$(FCMODEXT)
endif
#
mo_cdi.f90: $(top_srcdir)/src/cdi.h
	ruby $(top_srcdir)/interfaces/f2003/bindGen.rb $(top_srcdir)/src/cdi.h $@
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
#CLEANFILES += cdilib.c
if CREATE_ISOC
CLEANFILES += $(top_builddir)/src/mo_cdi.$(FCMODEXT) $(top_builddir)/src/mo_cdi.o
endif

if ENABLE_CDI_LIB
CLEANFILES += pkgconfig/cdi.pc

install-exec-local: pkgconfig/cdi.pc
	$(mkinstalldirs) "$(DESTDIR)$(libdir)/pkgconfig"
	$(install_sh_DATA) pkgconfig/cdi.pc "$(DESTDIR)$(libdir)/pkgconfig/cdi.pc"

uninstall-local:
	rm -f "$(DESTDIR)$(libdir)/pkgconfig/cdi.pc"
	rmdir "$(DESTDIR)$(libdir)/pkgconfig"

endif

install-exec-hook:
	-@rmdir "$(DESTDIR)$(libdir)"
install-data-hook:
	-@rmdir "$(DESTDIR)$(includedir)"

