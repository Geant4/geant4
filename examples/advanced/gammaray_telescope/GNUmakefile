# $Id: GNUmakefile,v 1.21 2002-12-08 18:17:25 flongo Exp $
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

ifdef G4STORE_DATA
  CPPFLAGS += -DG4STORE_DATA
endif

ifdef G4ANALYSIS_USE
  CPPFLAGS += `aida-config --include`
  LDFLAGS += `aida-config --lib`
endif

ifdef G4ANALYSIS_USE_NTUPLE
  CPPFLAGS += -DG4ANALYSIS_USE_NTUPLE
endif

.PHONY: all 
all:  lib bin

include $(G4INSTALL)/config/binmake.gmk

