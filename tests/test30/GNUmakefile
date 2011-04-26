# $Id: GNUmakefile,v 1.13 2010-01-12 09:28:48 gcosmo Exp $
# --------------------------------------------------------------
# GNUmakefile for tests module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := test30
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
