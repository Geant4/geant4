# $Id: GNUmakefile,v 1.6 2001-12-07 14:18:01 pfeiffer Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := XrayFluo
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../../
endif

# Copy the GammaRayTelAnalysis files according to which version
# is selected. For this the target tmpana has been created.
# As the dependency is checked before "all" is reached, we need
# to do some tricks here ...

ifdef G4ANALYSIS_AIDA_CONFIG_CFLAGS
# AIDA (version 2.2 for now) selected
FOO := $(shell cp include/XrayFluoAnalysisManager.hh-AIDA include/XrayFluoAnalysisManager.hh )
FOO := $(shell cp src/XrayFluoAnalysisManager.cc-AIDA src/XrayFluoAnalysisManager.cc )
else
# Anaphe selected 
include anaphe.gmk
FOO := $(shell cp include/XrayFluoAnalysisManager.hh-Anaphe include/XrayFluoAnalysisManager.hh )
FOO := $(shell cp src/XrayFluoAnalysisManager.cc-Anaphe src/XrayFluoAnalysisManager.cc )
endif

# End Analysis.


.PHONY: all
all: tmpana lib bin 

include $(G4INSTALL)/config/architecture.gmk

# include the configuration for the Analysis
include $(G4INSTALL)/config/analysis.gmk

include $(G4INSTALL)/config/binmake.gmk

visclean:
	rm -f g4*.prim g4*.eps g4*.wrl
	rm -f .DAWN_*

