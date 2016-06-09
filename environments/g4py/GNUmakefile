# $Id: GNUmakefile,v 1.1 2006/02/27 09:29:42 kmura Exp $
# $Name: geant4-09-00 $
# ===========================================================
#   Makefile for building Geant4Py
# ===========================================================
include  ./config/config.gmk

SUBDIR := source site-modules

.PHONY: all install clean uninstall

all:
	@for dir in $(SUBDIR); do (cd $$dir && $(MAKE)); done;:

install:
	@for dir in $(SUBDIR); do (cd $$dir && $(MAKE) install); done;:

clean:
	@for dir in $(SUBDIR); do (cd $$dir && $(MAKE) clean); done;:

uninstall:
	@-\rm -f -r $(Q_LIBDIR)
	@-\rm -f -r $(Q_INCDIR)

