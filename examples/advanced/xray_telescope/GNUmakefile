# --------------------------------------------------------------
# $Id: GNUmakefile,v 1.7 2002-11-14 15:56:14 santin Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := XrayTel
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

G4VIS_USE = 1

include $(G4INSTALL)/config/architecture.gmk

ifdef G4ANALYSIS_USE_PLOTTER
  CPPFLAGS += -DG4ANALYSIS_USE_PLOTTER
endif

ifdef G4ANALYSIS_USE
  CPPFLAGS += `aida-config --include`
  LDFLAGS += `aida-config --lib`
endif

.PHONY: all
all: lib bin

include $(G4INSTALL)/config/binmake.gmk

