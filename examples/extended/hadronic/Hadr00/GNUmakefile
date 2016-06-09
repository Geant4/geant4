# $Id: GNUmakefile,v 1.2 2010-01-11 16:40:31 gcosmo Exp $
# GEANT4 tag $Name: not supported by cvs2svn $
# --------------------------------------------------------------
# GNUmakefile for examples module 
# --------------------------------------------------------------

name := Hadr00
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../../..
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

