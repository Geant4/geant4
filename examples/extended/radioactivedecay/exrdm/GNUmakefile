# $Id: GNUmakefile,v 1.6 2008-08-10 13:21:09 vnivanch Exp $
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


.PHONY: all
all: lib bin

include $(G4INSTALL)/config/binmake.gmk

