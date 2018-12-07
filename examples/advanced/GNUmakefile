# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

ifndef G4INSTALL
  G4INSTALL = ../..
endif
 
include $(G4INSTALL)/config/architecture.gmk

SUBDIRS =  gammaray_telescope underground_physics xray_fluorescence xray_telescope
SUBDIRS += brachytherapy hadrontherapy iort_therapy medical_linac
SUBDIRS += purging_magnet radioprotection human_phantom
SUBDIRS += air_shower microbeam microelectronics nanobeam
SUBDIRS += composite_calorimeter lAr_calorimeter
SUBDIRS += amsEcal ChargeExchangeMC eRosita

.PHONY : all clean clean_libs

all:
	@for dir in $(SUBDIRS); do (cd $$dir; $(MAKE)); done

clean:
	@for dir in $(SUBDIRS); do (cd $$dir; $(MAKE) clean); done

clean_libs:
	@for dir in $(SUBDIRS); do (cd $$dir; $(MAKE) clean_libs); done
