# $Id: GNUmakefile,v 1.1 2010-09-16 16:37:12 gcosmo Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module
# --------------------------------------------------------------

name := rdecay02
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../../..
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

anaclean:
	rm -f $(G4WORKDIR)/tmp/$(G4SYSTEM)/$(G4TARGET)/exrdmHisto* 
	rm -f $(G4WORKDIR)/tmp/$(G4SYSTEM)/$(G4TARGET)/exrdmAnalysis* 
	rm -f $(G4WORKDIR)/tmp/$(G4SYSTEM)/$(G4TARGET)/exrdmStepping* 
