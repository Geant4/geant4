# --------------------------------------------------------------
# $Id: GNUmakefile,v 1.4 2002-06-10 16:09:04 guatelli Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := Brachy

G4TARGET := $(name)
G4EXLIB := true


G4ANALYSIS_USE = 1


ifndef G4INSTALL
  G4INSTALL = ../../..
endif

CPPFLAGS += -I${LHCXX_REL_DIR}/include -I/afs/cern.ch/sw/contrib/AIDA/2.2/cpp

.PHONY: all
all: lib bin


include $(G4INSTALL)/config/binmake.gmk


