# $Id: GNUmakefile,v 1.1 2003-10-08 16:31:45 hpw Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := exrdm01
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../../../../../../..
endif

.PHONY: all
all: lib bin

CPPFLAGS +=-I$(G4INSTALL)/source/processes/hadronic/models/radiative_decay/include
EXTRALIBS +=-lG4hadronic_radioactivedecay -lG4hadronic_deex -lG4hadronic_util_gen

include $(G4INSTALL)/config/binmake.gmk

