# --------------------------------------------------------------
# $Id: GNUmakefile,v 1.3 2002-12-16 17:36:39 guatelli Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := test50

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


