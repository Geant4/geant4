# $Id: GNUmakefile,v 1.14 2001-12-07 16:48:53 gcosmo Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := GammaRayTel
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

include $(G4INSTALL)/config/architecture.gmk

# Copy the GammaRayTelAnalysis files according to which version
# is selected.
# As the dependency is checked before "all" is reached, we need
# to do some tricks here ...

ifdef G4ANALYSIS_AIDA_CONFIG_CFLAGS
# AIDA (version 2.2 for now) selected
  FOO := $(shell cp include/GammaRayTelAnalysis.hh-AIDA include/GammaRayTelAnalysis.hh )
  FOO := $(shell cp src/GammaRayTelAnalysis.cc-AIDA src/GammaRayTelAnalysis.cc )
else
# Anaphe selected 
  include anaphe.gmk
  FOO := $(shell cp include/GammaRayTelAnalysis.hh-Anaphe include/GammaRayTelAnalysis.hh )
  FOO := $(shell cp src/GammaRayTelAnalysis.cc-Anaphe src/GammaRayTelAnalysis.cc )
endif

# End Analysis.

.PHONY: all 
all:  lib bin

include $(G4INSTALL)/config/binmake.gmk

