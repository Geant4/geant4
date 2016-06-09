# $Id: GNUmakefile,v 1.12 2009/08/13 20:48:04 kaitanie Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := Hadrontherapy
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

.PHONY: all
all: lib bin

include $(G4INSTALL)/config/architecture.gmk

ifdef G4ANALYSIS_USE
CPPFLAGS += -DANALYSIS_USE
endif
ifndef G4ANALYSIS_USE      # If we don't have AIDA
ifdef G4ANALYSIS_USE_ROOT   # And we have ROOT
CPPFLAGS += -DANALYSIS_USE -DG4ANALYSIS_USE_ROOT
CPPFLAGS += $(shell root-config --cflags)
LDFLAGS  += $(shell root-config --glibs)
endif
endif

include $(G4INSTALL)/config/binmake.gmk

ifdef G4ANALYSIS_USE 
endif
