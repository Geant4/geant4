# $Id: GNUmakefile,v 1.2 2006/04/25 10:17:30 kmura Exp $
# $Name: geant4-08-02 $
# ===========================================================
#   Makefile for building Geant4Py site-modules
# ===========================================================

#SUBDIR := EZdetector EZgeom python
SUBDIR  := python
SUBDIR2 := EZgeom

.PHONY: all install clean cleanlib

all:
	@for dir in $(SUBDIR); do (cd $$dir && $(MAKE)); done;:
	@for dir in $(SUBDIR2); do (cd $$dir && $(MAKE)); done;:

install:
	@for dir in $(SUBDIR); do (cd $$dir && $(MAKE) install); done;:
	@for dir in $(SUBDIR2); do (cd $$dir && $(MAKE) installall); done;:

clean:
	@for dir in $(SUBDIR); do (cd $$dir && $(MAKE) clean); done;:
	@for dir in $(SUBDIR2); do (cd $$dir && $(MAKE) cleanall); done;:

