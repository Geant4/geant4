# --------------------------------------------------------------
# $Id: GNUmakefile,v 1.9 2010-01-11 16:20:08 gcosmo Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := XrayTel
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../..
endif

include $(G4INSTALL)/config/architecture.gmk

ifdef G4ANALYSIS_USE_PLOTTER
  CPPFLAGS += -DG4ANALYSIS_USE_PLOTTER
endif

.PHONY: all
all: lib bin

include $(G4INSTALL)/config/binmake.gmk

