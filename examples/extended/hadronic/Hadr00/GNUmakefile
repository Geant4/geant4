# $Id: GNUmakefile,v 1.1 2008/07/07 16:37:26 vnivanch Exp $
# GEANT4 tag $Name: geant4-09-02 $
# --------------------------------------------------------------
# GNUmakefile for examples module 
# --------------------------------------------------------------

name := hadr00
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../../
endif

.PHONY: all
all: lib bin

include $(G4INSTALL)/config/architecture.gmk

include $(G4INSTALL)/config/binmake.gmk

histclean:
	rm -f $(G4WORKDIR)/tmp/$(G4SYSTEM)/$(G4TARGET)/Histo.o

visclean:
	rm -f g4*.prim g4*.eps g4*.wrl
	rm -f .DAWN_*

