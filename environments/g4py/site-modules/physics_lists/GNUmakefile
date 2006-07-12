# $Id: GNUmakefile,v 1.3 2006-07-12 08:45:14 kmura Exp $
# $Name: not supported by cvs2svn $
# ===========================================================
#   Makefile for building Geant4Py site-modules
# ===========================================================

SUBDIR := ExN01pl ExN03pl EMSTDpl EMLowEpl
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

