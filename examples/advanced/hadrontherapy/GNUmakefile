# $Id: GNUmakefile 68008 2013-03-13 11:30:30Z gcosmo $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := hadrontherapy
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
EXTRALIBS  += $(shell root-config --glibs)
endif

include $(G4INSTALL)/config/binmake.gmk

