# $Id: GNUmakefile,v 1.1 2007-05-24 21:57:02 tinslay Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := BremsstrahlungSplitting_G4WrapperProcess
G4TARGET := $(name)
G4EXLIB := true
G4WORKDIR := .

ifndef G4INSTALL
  G4INSTALL = ../../../../../../../
endif

.PHONY: all
all: lib bin

include $(G4INSTALL)/config/binmake.gmk
