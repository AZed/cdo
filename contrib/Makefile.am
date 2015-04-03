EXTRA_DIST = makecompl.rb cdoCompletion.bash cdoCompletion.tcsh cdoCompletion.zsh

all-local:
	if test -x /usr/bin/ruby -o -x /usr/local/bin/ruby ; then \
	  ruby $(srcdir)/makecompl.rb -o completions4cdo; \
	fi

CLEANFILES  = `ls completions4cdo.*`
