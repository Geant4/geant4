# --------------------------------------------------------------
# $Id: GNUmakefile,v 1.9 2003-03-12 17:21:21 pia Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := test50
G4TARGET := $(name)
G4EXLIB := true

#G4ANALYSIS_USE = 1

ifndef G4INSTALL
  G4INSTALL = ../..
endif

include $(G4INSTALL)/config/architecture.gmk

CPPFLAGS += `aida-config --include`
LDFLAGS += `aida-config --lib`

#CPPFLAGS += -I${LHCXX_REL_DIR}/include -I/afs/cern.ch/sw/contrib/AIDA/3.0/cpp -g

.PHONY: all
all: lib bin

include $(G4INSTALL)/config/binmake.gmk


