# $Id: GNUmakefile,v 1.4 2002-08-29 15:37:24 dressel Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := exampleB01
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

.PHONY: all
all: lib bin

include /afs/cern.ch/user/d/dressel/G4scoring/scoring.gmk

CPPFLAGS +=

include $(G4INSTALL)/config/binmake.gmk
