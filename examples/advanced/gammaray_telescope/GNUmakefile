# $Id: GNUmakefile,v 1.23 2006/06/02 05:24:20 flongo Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := GammaRayTel
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../..
endif
.PHONY: all 
all:  lib bin

include $(G4INSTALL)/config/binmake.gmk

include hadronic_lists.gmk

ifdef G4ANALYSIS_USE
  CPPFLAGS += `aida-config --include`
  LDFLAGS += `aida-config --lib`
  LOADLIBS += `aida-config --lib`
endif

ifdef G4STORE_DATA
  CPPFLAGS += -DG4STORE_DATA
endif
