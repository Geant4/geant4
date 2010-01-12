# $Id: GNUmakefile,v 1.4 2010-01-12 09:43:13 gcosmo Exp $
# --------------------------------------------------------------
# GNUmakefile for tests module with AIDA interface
# --------------------------------------------------------------

name := test45
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../..
endif

.PHONY: all
all: lib bin

ifdef G4ANALYSIS_USE_ROOT
  CPPFLAGS += -DG4ANALYSIS_USE_ROOT
endif

include $(G4INSTALL)/config/architecture.gmk

ifdef G4ANALYSIS_USE_ROOT
  CPPFLAGS += $(shell $(ROOTSYS)/bin/root-config --cflags)
  LDFLAGS  += $(shell $(ROOTSYS)/bin/root-config --glibs)
endif

include $(G4INSTALL)/config/binmake.gmk

histclean:
	rm -f $(G4WORKDIR)/tmp/$(G4SYSTEM)/$(G4TARGET)/Histo.o

visclean:
	rm -f g4*.prim g4*.eps g4*.wrl
	rm -f .DAWN_*

