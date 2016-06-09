# $Id: GNUmakefile,v 1.10 2007/12/05 08:25:56 gcosmo Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.
# --------------------------------------------------------------

name := load_gdml
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../..
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
