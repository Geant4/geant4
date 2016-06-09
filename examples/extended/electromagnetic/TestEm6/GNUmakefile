# $Id: GNUmakefile,v 1.12 2004/09/27 09:35:30 maire Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := TestEm6
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../../..
endif

.PHONY: all
all: lib bin

#### G4ANALYSIS_USE := true

ifdef G4ANALYSIS_USE
  CPPFLAGS += -DG4ANALYSIS_USE
endif

include $(G4INSTALL)/config/architecture.gmk

ifdef G4ANALYSIS_USE
  # for the aida-config command see the README file
  CPPFLAGS += `aida-config --include`
  LDFLAGS  += `aida-config --lib`
endif

include $(G4INSTALL)/config/binmake.gmk

histclean:
	rm -f $(G4WORKDIR)/tmp/$(G4SYSTEM)/$(G4TARGET)/RunAction.o
	rm -f $(G4WORKDIR)/tmp/$(G4SYSTEM)/$(G4TARGET)/SteppingAction.o

visclean:
	rm -f g4*.prim g4*.eps g4*.wrl
	rm -f .DAWN_*
