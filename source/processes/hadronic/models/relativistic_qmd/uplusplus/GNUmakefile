# $Id: GNUmakefile,v 1.1 2000-12-12 16:57:55 hweber Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := Uplusplus
G4TARGET := $(name)
G4EXLIB = true
G4DEBUG = 1

ifndef G4INSTALL
  G4INSTALL = ../../../../../../../
endif

.PHONY: all
all: lib bin

# G4EXEC_BUILD = true
CPPFLAGS +=-I$(G4INSTALL)/source/processes/hadronic/models/relativistic_qmd/include
EXTRALIBS +=-lG4hadronic_urqmd
# LDLIBS +=-lG4hadronic_urqmd

include $(G4INSTALL)/config/binmake.gmk

