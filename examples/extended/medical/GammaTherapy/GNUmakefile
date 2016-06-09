# $Id: GNUmakefile,v 1.3 2008/06/11 15:17:57 vnivanch Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := GammaTherapy
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../geant4
endif

.PHONY: all
all: lib bin

include $(G4INSTALL)/config/architecture.gmk

include $(G4INSTALL)/config/binmake.gmk

histclean:
	rm -f $(G4WORKDIR)/tmp/$(G4SYSTEM)/$(G4TARGET)/Histo.o

visclean:
	rm -f $(G4WORKDIR)/tmp/$(G4SYSTEM)/$(G4TARGET)/VisManager.o
