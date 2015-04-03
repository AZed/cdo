EXTRA_DIST = makecompl.rb cdoCompletion.bash cdoCompletion.tcsh cdoCompletion.zsh
CLEANFILES = 

completions:
	if test -x /usr/bin/ruby -o -x /usr/local/bin/ruby ; then \
	  ruby -KN $(srcdir)/makecompl.rb -o cdoCompletion; \
	fi

if MAINTAINER_MODE
all-local: completions

CLEANFILES  += `ls cdoCompletion.*`
endif
