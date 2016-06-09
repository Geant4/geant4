# $Id: GNUmakefile,v 1.19 2004/06/23 09:33:55 stesting Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := TestEm2
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../../..
endif

.PHONY: all
all: lib bin

# G4ANALYSIS_USE := true

ifdef G4ANALYSIS_USE
  USE_AIDA := true
  CPPFLAGS += -DG4ANALYSIS_USE
endif

ifdef USE_AIDA
  CPPFLAGS += -DUSE_AIDA
endif

include $(G4INSTALL)/config/architecture.gmk

ifdef USE_AIDA
  # for the aida-config command see the README file
  CPPFLAGS += `aida-config --include`
  LDFLAGS  += `aida-config --lib`
endif

include $(G4INSTALL)/config/binmake.gmk

histclean:
	rm -f $(G4WORKDIR)/tmp/$(G4SYSTEM)/$(G4TARGET)/RunAction.o

visclean:
	rm -f g4*.prim g4*.eps g4*.wrl
	rm -f .DAWN_*
