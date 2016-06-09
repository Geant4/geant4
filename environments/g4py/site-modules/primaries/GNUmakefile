# $Id: GNUmakefile,v 1.2 2006/04/25 10:29:52 kmura Exp $
# $Name: geant4-09-01 $
# ===========================================================
#   Makefile for building Geant4Py site-modules
# ===========================================================

SUBDIR := ParticleGun MedicalBeam

.PHONY: all install clean cleanlib

all:
	@for dir in $(SUBDIR); do (cd $$dir && $(MAKE)); done;:

install:
	@for dir in $(SUBDIR); do (cd $$dir && $(MAKE) install); done;:

clean:
	@for dir in $(SUBDIR); do (cd $$dir && $(MAKE) clean); done;:

