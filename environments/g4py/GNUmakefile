# $Id$
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

