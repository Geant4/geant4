# --------------------------------------------------------------
# $Id: GNUmakefile,v 1.8 2003-05-09 07:19:30 gcosmo Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := Brachy

G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL =~/GEANT4/geant4/
endif

ifdef G4ANALYSIS_USE  
  CPPFLAGS += -I${LHCXX_REL_DIR}/include -I/afs/cern.ch/sw/contrib/AIDA/3.0/cpp -g
endif

.PHONY: all
all: lib bin


include $(G4INSTALL)/config/binmake.gmk


