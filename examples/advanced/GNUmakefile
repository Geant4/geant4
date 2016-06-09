# $Id: GNUmakefile,v 1.5 2004/06/23 16:09:48 gcosmo Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

ifndef G4INSTALL
  G4INSTALL = ../..
endif
 
include $(G4INSTALL)/config/architecture.gmk

SUBDIRS =  gammaray_telescope underground_physics xray_fluorescence xray_telescope
SUBDIRS += brachytherapy cosmicray_charging hadrontherapy
SUBDURS += medical_linac purging_magnet radioprotection
SUBDIRS += Tiara Rich composite_calorimeter lAr_calorimeter

.PHONY : all clean clean_libs

all:
	@for dir in $(SUBDIRS); do (cd $$dir; $(MAKE)); done

clean:
	@for dir in $(SUBDIRS); do (cd $$dir; $(MAKE) clean); done

clean_libs:
	@for dir in $(SUBDIRS); do (cd $$dir; $(MAKE) clean_libs); done
