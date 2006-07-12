# $Id: GNUmakefile,v 1.3 2006-07-12 01:56:23 kmura Exp $
# $Name: not supported by cvs2svn $
# ===========================================================
#   Makefile for building Geant4Py site-modules
# ===========================================================
include  ../config/config.gmk

SUBDIR := geometries materials physics_lists primaries utils
SUBDIR += EZsim

.PHONY: all install clean uninstall

all:
	@for dir in $(SUBDIR); do (cd $$dir && $(MAKE)); done;:

install:
	@for dir in $(SUBDIR); do (cd $$dir && $(MAKE) install); done;:

clean:
	@for dir in $(SUBDIR); do (cd $$dir && $(MAKE) clean); done;:

uninstall:
	@-\rm -f -r $(Q_LIBDIR)/site-modules
	@-\rm -f $(Q_LIBDIR)/lib*.so


