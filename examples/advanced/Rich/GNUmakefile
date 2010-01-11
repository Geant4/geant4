# $Id: GNUmakefile,v 1.3 2010-01-11 14:24:55 gcosmo Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := RichTbSim
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../..
endif

ifdef G4ANALYSIS_USE
  CPPFLAGS += -DG4ANALYSIS_USE
endif


include $(G4INSTALL)/config/architecture.gmk

#
.PHONY: all
all: lib bin

include $(G4INSTALL)/config/binmake.gmk

visclean:
	rm -f g4*.prim g4*.eps g4*.wrl
	rm -f .DAWN_*












