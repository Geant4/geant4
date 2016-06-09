# $Id: GNUmakefile,v 1.6 2003/11/13 12:00:27 vnivanch Exp $
# GEANT4 tag $Name: geant4-06-00 $
# --------------------------------------------------------------
# GNUmakefile for examples module with AIDA interface
# --------------------------------------------------------------

name := TestEm9
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../../../
endif

.PHONY: all
all: lib bin

ifdef G4ANALYSIS_USE
  CPPFLAGS += -DG4ANALYSIS_USE
endif

include $(G4INSTALL)/config/architecture.gmk

# for the aida-config command see the README file
ifdef G4ANALYSIS_USE
  CPPFLAGS += `aida-config --include`
  LDFLAGS  += `aida-config --lib`
endif

include $(G4INSTALL)/config/binmake.gmk

analclean:
	rm -f $(G4WORKDIR)/tmp/$(G4SYSTEM)/$(G4TARGET)/Histo.o

visclean:
	rm -f $(G4WORKDIR)/tmp/$(G4SYSTEM)/$(G4TARGET)/VisManager.o


