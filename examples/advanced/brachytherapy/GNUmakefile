# --------------------------------------------------------------
# $Id: GNUmakefile,v 1.7 2002-11-15 17:36:55 guatelli Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := Brachy

G4TARGET := $(name)
G4EXLIB := true


G4ANALYSIS_USE = 1


ifndef G4INSTALL
  G4INSTALL =~/GEANT4/geant4/
endif

CPPFLAGS += -I${LHCXX_REL_DIR}/include -I/afs/cern.ch/sw/contrib/AIDA/3.0/cpp -g


.PHONY: all
all: lib bin


include $(G4INSTALL)/config/binmake.gmk


