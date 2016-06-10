# $Id: GNUmakefile 66508 2012-12-19 10:16:45Z gcosmo $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := GammaRayTel
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../..
endif

.PHONY: all
all:  had_list lib bin

had_list:
	@echo "Gamma-Ray Telescope advanced example"

include $(G4INSTALL)/config/binmake.gmk

