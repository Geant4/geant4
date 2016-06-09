# $Id: GNUmakefile,v 1.8 2008/12/03 07:33:39 kmura Exp $
# $Name: geant4-09-02 $
# ===========================================================
#   Makefile for building Geant4Py modules
# ===========================================================
include  ../config/config.gmk

SUBDIR1 := global interface intercoms run event tracking track
SUBDIR1 += geometry materials particles processes digits_hits
SUBDIR2 := visualization graphics_reps physics_lists gdml
SUBDIR3 := python

.PHONY: all install clean 

all:
	@for dir in $(SUBDIR1); do (cd $$dir && $(MAKE)); done;:
	@for dir in $(SUBDIR2); do (cd $$dir && $(MAKE)); done;:
	@for dir in $(SUBDIR3); do (cd $$dir && $(MAKE)); done;:

install:
	@for dir in $(SUBDIR1); do (cd $$dir && $(MAKE) install); done;:
	@for dir in $(SUBDIR2); do (cd $$dir && $(MAKE) install); done;:
	@for dir in $(SUBDIR3); do (cd $$dir && $(MAKE) install); done;:

clean:
	@for dir in $(SUBDIR1); do (cd $$dir && $(MAKE) clean); done;:
	@for dir in $(SUBDIR2); do (cd $$dir && $(MAKE) clean); done;:
	@for dir in $(SUBDIR3); do (cd $$dir && $(MAKE) clean); done;:

