# $Id: GNUmakefile,v 1.2 2000-08-08 12:27:09 hweber Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := Uplusplus
G4TARGET := $(name)

G4DEBUG = 1

ifndef G4INSTALL
  G4INSTALL = ../../../../../../../
endif

G4EXEC_BUILD = true
CPPFLAGS +=-I$(G4INSTALL)/source/processes/hadronic/models/relativistic_qmd/include
EXTRALIBS +=-lG4hadronic_urqmd
# LDLIBS +=-lG4hadronic_urqmd

include $(G4INSTALL)/config/binmake.gmk

