# $Id: GNUmakefile,v 1.2 2006/07/05 11:41:20 gcosmo Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------
name := hepmcEx02-clhep
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../../../..
endif

.PHONY: all
all : lib bin

include $(G4INSTALL)/config/binmake.gmk

# Path for PYTHIA Fortran library. Based on CERNLIB-2005.
# Add "/cern/pro/bin" to $PATH first !
#
LDLIBS1     += $(shell cernlib -v pro pythia6205 pdflib804 mathlib) -lg2c
