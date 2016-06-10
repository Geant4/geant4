# $Id: GNUmakefile 76296 2013-11-08 13:17:47Z gcosmo $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := ParN04
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../../..
endif

include ParGNUmakefile

.PHONY: all
all: lib bin

include $(G4INSTALL)/config/binmake.gmk
