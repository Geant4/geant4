# $Id: GNUmakefile,v 1.2 2008/08/27 11:04:02 gcosmo Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.
# --------------------------------------------------------------

name := load_gdml
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../../../..
endif

.PHONY: all
all: lib bin
ifndef XERCESCROOT	
	@echo XERCESCROOT not defined!
endif

include $(G4INSTALL)/config/binmake.gmk

visclean:
	rm -f g4*.prim g4*.eps g4*.wrl
	rm -f .DAWN_*
