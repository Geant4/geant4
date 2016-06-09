# $Id: GNUmakefile,v 1.2 2006/07/05 13:07:54 gcosmo Exp $
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

include hadronic_lists.gmk

.PHONY: all
all: had_list lib bin

had_list:
	@echo "This example uses the QGSP 'educated guess physics list'"
	@echo "(see $(G4INSTALL)/physics_lists/hadronic)."
	@echo "The hadronic physics lists should be compiled before"
	@echo "compiling this example."

include $(G4INSTALL)/config/binmake.gmk
