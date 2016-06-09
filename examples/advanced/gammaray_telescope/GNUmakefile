# $Id: GNUmakefile,v 1.24 2006/11/17 15:49:35 flongo Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := GammaRayTel
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

.PHONY: all
all:  had_list lib bin

had_list:
	@echo "Gamma Ray Tel Advanced example"

include $(G4INSTALL)/config/binmake.gmk

