# --------------------------------------------------------------
# $Id: GNUmakefile,v 1.5 2009-09-15 16:34:42 cirrone Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := MedLinac

G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

.PHONY: all
all: lib bin
include $(G4INSTALL)/config/binmake.gmk


ifdef G4STORE_DATA
  CPPFLAGS += -DG4STORE_DATA
endif

