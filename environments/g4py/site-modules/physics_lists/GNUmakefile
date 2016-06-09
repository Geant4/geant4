# $Id: GNUmakefile,v 1.2 2006/04/25 10:31:40 kmura Exp $
# $Name: geant4-08-01 $
# ===========================================================
#   Makefile for building Geant4Py site-modules
# ===========================================================

SUBDIR := ExN01pl ExN03pl EMSTDpl 
#SUBDIR2 := GenericPhysicsList

.PHONY: all install clean

all:
	@for dir in $(SUBDIR); do (cd $$dir && $(MAKE)); done;:
	#@for dir in $(SUBDIR2); do (cd $$dir && $(MAKE)); done;:

install:
	@for dir in $(SUBDIR); do (cd $$dir && $(MAKE) install); done;:
	#@for dir in $(SUBDIR2); do (cd $$dir && $(MAKE) installall); done;:

clean:
	@for dir in $(SUBDIR); do (cd $$dir && $(MAKE) clean); done;:
	#@for dir in $(SUBDIR2); do (cd $$dir && $(MAKE) cleanall); done;:

