# $Id: GNUmakefile,v 1.1 2002-12-13 15:47:31 stesting Exp $
# --------------------------------------------------------------
# GNUmakefile for tests module.  Gabriele Cosmo, 27/06/98.
# --------------------------------------------------------------

name := Marker
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../..
endif

.PHONY: all

all: lib bin

include $(G4INSTALL)/config/scripts/binmake.gmk

CXXFLAGS_FINAL := $(CXXFLAGS) 
INCFLAGS_FINAL := $(INCFLAGS)
LDLIBS_FINAL := $(LDLIBS)
#LDFLAGS_FINAL := $(LDFLAGS)
