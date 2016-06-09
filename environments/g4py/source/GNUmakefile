# $Id: GNUmakefile,v 1.2 2006/04/25 08:09:43 kmura Exp $
# $Name: geant4-08-01 $
# ===========================================================
#   Makefile for building Geant4Py modules
# ===========================================================
include  ../config/config.gmk

SUBDIR1 := global interface run event tracking track
SUBDIR1 += geometry materials particles processes digits_hits
SUBDIR2 := visualization graphics_reps 
SUBDIR3 := python
SUBDIR4 := boost version

.PHONY: all install clean uninstall

all:
	@for dir in $(SUBDIR1); do (cd $$dir && $(MAKE)); done;:
	@for dir in $(SUBDIR2); do (cd $$dir && $(MAKE)); done;:
	@for dir in $(SUBDIR3); do (cd $$dir && $(MAKE)); done;:

install:
	@for dir in $(SUBDIR1); do (cd $$dir && $(MAKE) install); done;:
	@for dir in $(SUBDIR2); do (cd $$dir && $(MAKE) install); done;:
	@for dir in $(SUBDIR3); do (cd $$dir && $(MAKE) install); done;:
	@for dir in $(SUBDIR4); do (cd $$dir && $(MAKE) install); done;:

clean:
	@for dir in $(SUBDIR1); do (cd $$dir && $(MAKE) clean); done;:
	@for dir in $(SUBDIR2); do (cd $$dir && $(MAKE) clean); done;:
	@for dir in $(SUBDIR3); do (cd $$dir && $(MAKE) clean); done;:

uninstall:
	@-\rm -f -r $(Q_LIBDIR)/Geant4
	@-\rm -f -r $(Q_INCDIR)

