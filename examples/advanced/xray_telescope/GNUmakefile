# --------------------------------------------------------------
# $Id: GNUmakefile,v 1.1 2000-10-04 16:38:30 nartallo Exp $
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
G4VIS_USE_DAWNFILE = 1	
G4VIS_USE_DAWN = 1	
G4VIS_USE_VRMLFILE = 1	
G4VIS_USE_VRML = 1

.PHONY: all
all: lib bin

include $(G4INSTALL)/config/binmake.gmk
