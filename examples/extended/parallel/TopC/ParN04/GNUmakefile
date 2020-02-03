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
