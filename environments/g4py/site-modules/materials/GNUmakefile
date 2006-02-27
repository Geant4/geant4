# $Id: GNUmakefile,v 1.1 2006-02-27 09:48:56 kmura Exp $
# $Name: not supported by cvs2svn $
# ===========================================================
#   Makefile for building Geant4Py site-modules
# ===========================================================

SUBDIR := Qmaterials NISTmaterials

.PHONY: all install clean cleanlib

all:
	@for dir in $(SUBDIR); do (cd $$dir && $(MAKE)); done;:

install:
	@for dir in $(SUBDIR); do (cd $$dir && $(MAKE) install); done;:

clean:
	@for dir in $(SUBDIR); do (cd $$dir && $(MAKE) clean); done;:

