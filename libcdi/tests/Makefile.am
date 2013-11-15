## Process this file with automake to produce Makefile.in
#
TESTS = test_grib.sh cksum_verify \
	test_cksum_grib test_cksum_nc test_cksum_extra \
	test_cksum_service test_cksum_nc2 test_cksum_nc4 test_cksum_ieg \
	test_chunk_cksum \
	pio_write_run pio_cksum_mpinonb pio_cksum_fpguard \
	pio_cksum_asynch pio_cksum_writer pio_cksum_cdf \
	test_resource_copy
check_PROGRAMS = cksum_verify test_grib cksum_write cksum_read pio_write \
	test_resource_copy cksum_write_chunk
#
test_grib_SOURCES = test_grib.c
cksum_verify_SOURCES = cksum_verify.c cksum.c cksum.h
cksum_write_SOURCES = cksum_write.c cksum.c cksum.h
cksum_write_chunk_SOURCES = cksum_write_chunk.c cksum.c cksum.h
cksum_read_SOURCES = cksum_read.c \
	var_cksum.c var_cksum.h \
	stream_cksum.c stream_cksum.h \
	cksum.c cksum.h \
	ensure_array_size.h ensure_array_size.c
pio_write_SOURCES = pio_write.c cksum.h cksum.c
test_resource_copy_SOURCES = test_resource_copy.c
#
AM_CFLAGS = $(YAXT_CFLAGS)
LDADD = ../src/libcdi.la -lm
INCLUDES   = -I$(top_srcdir)/src
#
EXTRA_DIST = $(TESTS)
#
CLEANFILES  = `ls *~ *.grb *.nc *.srv *.ext example_*.cksum`
#
AUTOMAKE_OPTIONS = color-tests
