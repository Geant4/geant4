# $Id: GNUmakefile,v 1.8 2004-12-13 16:11:53 gcosmo Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := exampleN04
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

include hadronic_lists.gmk

.PHONY: all
all:  had_list lib bin

had_list:
	@echo "This example uses the QGSP 'educated guess physics list'"
	@echo "(see $G4INSTALL/physics_lists/hadronic)."
	@echo "The hadronic physics lists should be compiled before"
	@echo "compiling this example."

include $(G4INSTALL)/config/binmake.gmk

