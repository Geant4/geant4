# $Id: GNUmakefile,v 1.3 2001-12-04 13:38:55 gcosmo Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

ifndef G4INSTALL
  G4INSTALL = ../..
endif
 
include $(G4INSTALL)/config/architecture.gmk

SUBDIRS =  brachytherapy gammaray_telescope underground_physics
SUBDIRS += xray_fluorescence xray_telescope

.PHONY : all clean clean_libs

all:
	@for dir in $(SUBDIRS); do (cd $$dir; $(MAKE)); done

clean:
	@for dir in $(SUBDIRS); do (cd $$dir; $(MAKE) clean); done

clean_libs:
	@for dir in $(SUBDIRS); do (cd $$dir; $(MAKE) clean_libs); done
