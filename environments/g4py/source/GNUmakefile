# $Id: GNUmakefile,v 1.9 2010-12-02 08:21:04 kmura Exp $
# $Name: not supported by cvs2svn $
# ===========================================================
#   Makefile for building Geant4Py modules
# ===========================================================
include  ../config/config.gmk

SUBDIR1 := global interface intercoms run event tracking track
SUBDIR1 += geometry materials particles processes digits_hits
SUBDIR2 := visualization graphics_reps physics_lists gdml

ifeq ($(Q_WITH_PYTHON3), 0)
  SUBDIR3 := python
else
  SUBDIR3 := python3
endif

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

