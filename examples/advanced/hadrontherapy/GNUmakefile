# $Id: GNUmakefile,v 1.16 2010-02-08 16:39:27 cirrone Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := Hadrontherapy
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../..
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
