# $Id: GNUmakefile,v 1.3 2006/11/16 13:38:46 gunter Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

include ParGNUmakefile

name := ParN04
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../..
endif


.PHONY: all
all: had_list lib bin

had_list:
	@echo "This example uses the QGSP 'educated guess physics list'"

include $(G4INSTALL)/config/binmake.gmk
