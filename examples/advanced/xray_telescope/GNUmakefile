# --------------------------------------------------------------
# $Id: GNUmakefile,v 1.4 2000-11-16 13:46:04 nartallo Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := XrayTel
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

G4VIS_USE = 1

.PHONY: all
all: lib bin

include $(G4INSTALL)/config/binmake.gmk

