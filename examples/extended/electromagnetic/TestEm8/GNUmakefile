# $Id: GNUmakefile,v 1.5 2003-11-24 18:09:22 vnivanch Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := TestEm8
G4TARGET := $(name)
G4EXLIB  := true

ifndef G4INSTALL
  G4INSTALL = ../../../..
endif

.PHONY: all
all: lib bin

include $(G4INSTALL)/config/architecture.gmk

include $(G4INSTALL)/config/binmake.gmk
