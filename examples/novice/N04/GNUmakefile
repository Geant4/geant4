# $Id: GNUmakefile,v 1.7 2004-12-02 14:28:03 santin Exp $
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
	@echo "This example (novice/N04) uses the QGSP 'educated guess physics list'"
	@echo "(see http://cmsdoc.cern.ch/~hpw/GHAD/HomePage)."
	@echo "The Hadronic Physics lists should be compiled before compiling this example." 
	@echo "If you haven't done so yet, please see the README file for detailed compilation instructions"

include $(G4INSTALL)/config/binmake.gmk

