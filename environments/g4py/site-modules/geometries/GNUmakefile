# $Id: GNUmakefile,v 1.3 2008/12/01 07:07:05 kmura Exp $
# $Name: geant4-09-02 $
# ===========================================================
#   Makefile for building Geant4Py site-modules
# ===========================================================

SUBDIR := Qgeom ExN01geom ExN03geom ezgeom

.PHONY: all install clean cleanlib

all:
	@for dir in $(SUBDIR); do (cd $$dir && $(MAKE)); done;:

install:
	@for dir in $(SUBDIR); do (cd $$dir && $(MAKE) install); done;:

clean:
	@for dir in $(SUBDIR); do (cd $$dir && $(MAKE) clean); done;:

