# $Id: GNUmakefile,v 1.9 2001-12-07 17:12:54 gcosmo Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := XrayFluo
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

include $(G4INSTALL)/config/architecture.gmk

# Copy the XrayFluoAnalysisManage files according to which version
# is selected.
# As the dependency is checked before "all" is reached, we need
# to do some tricks here ...

ifdef G4ANALYSIS_AIDA_CONFIG_CFLAGS
# AIDA (version 2.2 for now) selected
  FOO := $(shell cp -p include/XrayFluoAnalysisManager.hh-AIDA include/XrayFluoAnalysisManager.hh )
  FOO := $(shell cp -p src/XrayFluoAnalysisManager.cc-AIDA src/XrayFluoAnalysisManager.cc )
else
# Anaphe selected 
  include anaphe.gmk
  FOO := $(shell cp -p include/XrayFluoAnalysisManager.hh-Anaphe include/XrayFluoAnalysisManager.hh )
  FOO := $(shell cp -p src/XrayFluoAnalysisManager.cc-Anaphe src/XrayFluoAnalysisManager.cc )
endif

# End Analysis.


.PHONY: all visclean
all:  lib bin 

include $(G4INSTALL)/config/binmake.gmk

visclean:
	rm -f g4*.prim g4*.eps g4*.wrl
	rm -f .DAWN_*

