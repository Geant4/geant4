# --------------------------------------------------------------
# $Id: GNUmakefile,v 1.3 2002-10-09 09:54:05 sunanda Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := CompositeCalorimeter

G4TARGET := $(name)
G4EXLIB := true


G4ANALYSIS_USE = 1


ifndef G4INSTALL
  G4INSTALL = ../../../
endif

CPPFLAGS += -I${LHCXX_REL_DIR}/include -I/afs/cern.ch/sw/contrib/AIDA/2.2/cpp


.PHONY: all
all: lib bin

#EXTRALIBS = LHEP/lib/Linux-g++/libLHEP.so Packaging/lib/Linux-g++/libPackaging.a

include $(G4INSTALL)/config/binmake.gmk
# From lowenergy/test/GNUmakefile
CPPFLAGS += `aida-config --useHBook --include`
LDFLAGS  += `aida-config --useHBook --lib`
LOADLIBS  += `aida-config --useHBook --lib`


