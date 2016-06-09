# $Id: GNUmakefile,v 1.6.2.1 2008/04/23 09:26:58 gcosmo Exp $
# $Name: geant4-09-01-patch-03 $
# ===========================================================
#   Makefile for building Geant4Py modules
# ===========================================================
include  ../config/config.gmk

SUBDIR1 := global interface intercoms run event tracking track
SUBDIR1 += geometry materials particles processes digits_hits
SUBDIR2 := visualization graphics_reps physics_lists
SUBDIR3 := python
SUBDIR4 := boost

ifeq ($(Q_ENABLE_G4GDML), 1)
  SUBDIR2 += g4gdml
endif

ifeq ($(Q_ENABLE_GDML), 1)
  SUBDIR2 += gdml
endif

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

