# $Id: GNUmakefile,v 1.2 2006/07/05 09:50:45 gcosmo Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------
name := hepmcEx01-clhep
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../../../..
endif

include hadronic_lists.gmk

.PHONY: all
all : had_list lib bin

had_list:
	@echo "This example uses the QGSP 'educated guess physics list'"
	@echo "(see $(G4INSTALL)/physics_lists/hadronic)."
	@echo "The hadronic physics lists should be compiled before"
	@echo "compiling this example."

include $(G4INSTALL)/config/binmake.gmk

# Path for PYTHIA Fortran library. Based on CERNLIB-2005.
# Add "/cern/pro/bin" to $PATH first !
#
LDLIBS1     += $(shell cernlib -v pro pythia6205 pdflib804 mathlib) -lg2c
