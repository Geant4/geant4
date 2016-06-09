# $Id: GNUmakefile,v 1.7 2006/10/31 15:05:33 vnivanch Exp $
# GEANT4 tag $Name: geant4-09-01 $
# --------------------------------------------------------------
# GNUmakefile for examples module with AIDA interface
# --------------------------------------------------------------

name := hadr01
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../../
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
  LOADLIBS += `aida-config --lib`
endif

include $(G4INSTALL)/config/binmake.gmk

histclean:
	rm -f $(G4WORKDIR)/tmp/$(G4SYSTEM)/$(G4TARGET)/Histo.o

visclean:
	rm -f g4*.prim g4*.eps g4*.wrl
	rm -f .DAWN_*

