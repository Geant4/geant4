# $Id: GNUmakefile,v 1.1 2004-11-26 16:33:18 hpw Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := ExGflash
G4TARGET := $(name)
G4EXLIB := true


ifndef G4INSTALL
  G4INSTALL = ../../../..
endif

.PHONY: all
all: lib bin
include $(G4INSTALL)/config/architecture.gmk
include $(G4INSTALL)/config/binmake.gmk

