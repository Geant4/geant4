# $Id: GNUmakefile,v 1.18 2002-06-18 15:42:36 griccard Exp $
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

ifdef G4ANALYSIS_USE
  CPPFLAGS += `aida-config --useHBook --include`
  LDFLAGS += `aida-config --useHBook --lib`
endif

.PHONY: all 
all:  lib bin

include $(G4INSTALL)/config/binmake.gmk

