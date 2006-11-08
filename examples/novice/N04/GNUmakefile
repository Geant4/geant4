# $Id: GNUmakefile,v 1.10 2006-11-08 16:43:57 gunter Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := exampleN04
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

.PHONY: all
all:  had_list lib bin

had_list:
	@echo "This example uses the QGSP 'educated guess physics list'"

include $(G4INSTALL)/config/binmake.gmk

