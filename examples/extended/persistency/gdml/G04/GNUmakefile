# $Id: GNUmakefile 68025 2013-03-13 13:43:46Z gcosmo $
# --------------------------------------------------------------
# GNUmakefile for examples module.
# --------------------------------------------------------------

name := gdml_det
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
