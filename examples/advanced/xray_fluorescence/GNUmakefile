# $Id: GNUmakefile,v 1.5 2001-12-06 19:48:19 pfeiffer Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := XrayFluo
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../../
endif

.PHONY: all
all: tmpana lib bin 

include $(G4INSTALL)/config/architecture.gmk

# include the configuration for the Analysis
include $(G4INSTALL)/config/analysis.gmk


# Check whether to use AIDA (default) or not.
# Copy the XrayFluoAnalysisManager files according to which version
# is selected. For this the target tmpana has been created.
# This part of the makefile will be obsolete soon (see also below).

ifdef G4ANALYSIS_AIDA_CONFIG_CFLAGS
# AIDA (version 2.2 for now) selected
tmpana:
	cp include/XrayFluoAnalysisManager.hh-AIDA include/XrayFluoAnalysisManager.hh
	cp src/XrayFluoAnalysisManager.cc-AIDA src/XrayFluoAnalysisManager.cc
else
# Anaphe selected (will be obsolete as soon as Anaphe is AIDA 2.2 compliant)
include anaphe.gmk
tmpana:
	cp include/XrayFluoAnalysisManager.hh-Anaphe include/XrayFluoAnalysisManager.hh
	cp src/XrayFluoAnalysisManager.cc-Anaphe src/XrayFluoAnalysisManager.cc
endif

# End check for AIDA.


include $(G4INSTALL)/config/binmake.gmk

visclean:
	rm -f g4*.prim g4*.eps g4*.wrl
	rm -f .DAWN_*

