# --------------------------------------------------------------
# $Id: GNUmakefile,v 1.2 2009-03-06 09:22:48 gbarrand Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := cell
G4TARGET := $(name)
G4EXLIB := true

G4ANALYSIS_USE = 1

ifndef G4INSTALL
  G4INSTALL = ../..
endif

include $(G4INSTALL)/config/architecture.gmk

.PHONY: all
all: lib bin

include $(G4INSTALL)/config/binmake.gmk


