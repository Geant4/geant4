# $Id: GNUmakefile,v 1.2 2006/05/24 12:58:49 maire Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := TestEm16
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
  # add -v to see what is going on. 
  # add RunAction.o such that the linker finds AIDA_createAnalysisFactory
  # for the aida-config command see the README file
  CPPFLAGS += `aida-config --include`
  LDFLAGS  += $(G4WORKDIR)/tmp/$(G4SYSTEM)/TestEm16/RunAction.o `aida-config --lib`
endif

include $(G4INSTALL)/config/binmake.gmk

histclean:
	rm -f $(G4WORKDIR)/tmp/$(G4SYSTEM)/$(G4TARGET)/RunAction.o
	rm -f $(G4WORKDIR)/tmp/$(G4SYSTEM)/$(G4TARGET)/SteppingAction.o

visclean:
	rm -f g4*.prim g4*.eps g4*.wrl
	rm -f .DAWN_*
