# $Id: GNUmakefile,v 1.17 2010-05-18 11:14:39 cirrone Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := Hadrontherapy
G4TARGET := $(name)
G4EXLIB := true

# Debug info
#CPPFLAGS += -g 

ifndef G4INSTALL
  G4INSTALL = ../..
endif

.PHONY: all
all: lib bin

include $(G4INSTALL)/config/architecture.gmk

ifdef G4ANALYSIS_USE_ROOT   # If we have ROOT
CPPFLAGS += -DG4ANALYSIS_USE_ROOT
CPPFLAGS += $(shell root-config --cflags)
LDFLAGS  += $(shell root-config --glibs)
endif

include $(G4INSTALL)/config/binmake.gmk

