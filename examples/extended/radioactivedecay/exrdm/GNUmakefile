# $Id: GNUmakefile,v 1.5 2008-05-01 14:16:25 flei Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module
# --------------------------------------------------------------

name := exRDM
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../../..
endif

ifdef G4ANALYSIS_USE_ROOT
CPPFLAGS += -DG4ANALYSIS_USE_ROOT 
CPPFLAGS += -D_REENTRANT -I$(ROOTSYS)/include

EXTRALIBS += $(shell $(ROOTSYS)/bin/root-config --libs)

endif

ifdef G4ANALYSIS_USE_AIDA
   CPPFLAGS += -DG4ANALYSIS_USE_AIDA
   CPPFLAGS += `aida-config --include`
   LDLIBS += `aida-config --lib`
endif

.PHONY: all
all: lib bin

include $(G4INSTALL)/config/binmake.gmk

