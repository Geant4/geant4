# $Id: GNUmakefile,v 1.2 2006/04/25 10:18:58 kmura Exp $
# $Name: geant4-08-01 $
# ===========================================================
#   Makefile for building Geant4Py site-modules
# ===========================================================
include  ../config/config.gmk

SUBDIR := geometries materials physics_lists primaries
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


