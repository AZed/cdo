CLEANFILES =

export

TESTS = $(top_srcdir)/test/test_Gridarea.sh \
        $(top_srcdir)/test/test_Cat.sh \
        $(top_srcdir)/test/test_Spectral.sh \
        $(top_srcdir)/test/test_Timstat.sh \
        $(top_srcdir)/test/test_Select.sh \
        $(top_srcdir)/test/test_Vertint.sh

#        $(top_srcdir)/test/test_Remap.sh \
#       $(top_srcdir)/test/test_info.py
#	$(top_srcdir)/test/test_diff.py
#	$(top_srcdir)/test/test_Arith.py
#	$(top_srcdir)/test/test_Arithc.py
#	$(top_srcdir)/test/test_Selvar.py
#	$(top_srcdir)/test/test_intgridbil.py

EXTRA_DIST = $(TESTS) $(top_srcdir)/test/testStreams.py

CDO   = $(top_builddir)/src/cdo

DATAPATH   = $(top_srcdir)/test/data

PYTHONPATH = $(top_srcdir)/contrib/python:$(top_srcdir)/test

listTests:
	@echo $(TESTS)

test: $(TESTS)

.PHONY: $(TESTS)

CLEANFILES += `ls *.pyc`

AUTOMAKE_OPTIONS = color-tests
