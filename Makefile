# By default, this makefile only compiles for Linux.
# Before compiling for windows, run
#   
#   make distclean
#
# so that all object files get deleted and then rebuilt.
# Then run
#
#   make win
#
# I haven't managed to configure make to make separate 
# object files for windows and linux AND automatically detect
# dependencies.  The problem is that the automatic dependency
# detection relies on built-in rules.  
#
# See https://stackoverflow.com/a/2481326/3611932

RM=rm -f

all: nix 

nix: 
	make -f Makefile_nix64
win:
	make -f Makefile_win64

clean:
	$(RM) *.o .depend

distclean:
	$(RM) *.o *~ .depend
